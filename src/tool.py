"""DAQplugin-styled ChimeraX graphical front end for DaQRotFix."""

try:
    from importlib.metadata import PackageNotFoundError, version as package_version
except Exception:  # pragma: no cover - older bundled Python fallback
    PackageNotFoundError = Exception
    package_version = None

from chimerax.core.errors import UserError
from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from Qt.QtCore import Qt
from Qt.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .correction import residue_label
from .local_mdff import (
    DEFAULT_ANNEALING_UPDATES,
    DEFAULT_BACKBONE_SPRING_CONSTANT,
    DEFAULT_FINAL_MINIMIZATION_UPDATES,
    DEFAULT_FINAL_TEMPERATURE,
    DEFAULT_INITIAL_TEMPERATURE,
    DEFAULT_LOCAL_RADIUS,
    DEFAULT_OUTLIER_THRESHOLD,
    _outliers_in_residues,
    optimize_sidechains_with_local_mdff,
)


def _bundle_version():
    if package_version is None:
        return "0.6.1"
    try:
        return package_version("ChimeraX-DaQRotFix")
    except Exception:
        return "0.6.1"


def _adaptive_ui_scale():
    try:
        from Qt.QtGui import QGuiApplication

        screen = QGuiApplication.primaryScreen()
        height = screen.availableGeometry().height() if screen is not None else 0
    except Exception:
        height = 0
    if height <= 720:
        return 0.82
    if height <= 820:
        return 0.88
    if height <= 950:
        return 0.94
    return 1.0


_CAMERA_ATTRIBUTES = (
    "position",
    "field_of_view",
    "field_width",
    "eye_separation_scene",
    "eye_separation_pixels",
)


def _capture_view(session):
    """Capture camera-only state without snapshotting model coordinates."""
    view = session.main_view
    camera = view.camera
    camera_state = {}
    for attribute in _CAMERA_ATTRIBUTES:
        if hasattr(camera, attribute):
            value = getattr(camera, attribute)
            if attribute == "position":
                value = value.copy()
            camera_state[attribute] = value
    return {
        "camera": camera_state,
        "clip_planes": [plane.copy() for plane in view.clip_planes.planes()],
        "center_of_rotation": tuple(view.center_of_rotation),
        "center_of_rotation_method": view.center_of_rotation_method,
    }


def _restore_view(session, state):
    """Restore camera and clipping while leaving optimized coordinates intact."""
    view = session.main_view
    for attribute, value in state["camera"].items():
        setattr(view.camera, attribute, value)
    view.clip_planes.replace_planes(
        [plane.copy() for plane in state["clip_planes"]]
    )
    view.center_of_rotation = state["center_of_rotation"]
    view.center_of_rotation_method = state["center_of_rotation_method"]


class DaQRotFixTool(ToolInstance):
    SESSION_ENDURING = False
    SESSION_SAVE = False
    HELP = "help:command:correctrotamers"

    _HEADERS = (
        "Residue",
        "Selected\nRotamer",
        "Probability",
        "Local Outliers\n(Before → After)",
        "Result",
    )

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "DaQRotFix"
        self._running = False
        self._report = None
        self._threshold_initialized = False

        self.tool_window = MainToolWindow(self, close_destroys=True)
        self._build_ui()

        self._model_add_handler = session.triggers.add_handler(
            "add models", self._on_models_changed
        )
        self._model_remove_handler = session.triggers.add_handler(
            "remove models", self._on_models_changed
        )
        self.tool_window.manage(None)

    def delete(self):
        if hasattr(self, "_model_add_handler"):
            self.session.triggers.remove_handler(self._model_add_handler)
        if hasattr(self, "_model_remove_handler"):
            self.session.triggers.remove_handler(self._model_remove_handler)
        super().delete()

    def _build_ui(self):
        parent = self.tool_window.ui_area
        root = QWidget(parent)
        scale = _adaptive_ui_scale()

        def sp(value, minimum=0):
            return max(minimum, int(round(value * scale)))

        root.setStyleSheet(self._stylesheet(sp))
        outer = QVBoxLayout(root)
        outer.setContentsMargins(sp(6), sp(6), sp(6), sp(6))
        outer.setSpacing(sp(4))

        self.tabs = QTabWidget(root)
        outer.addWidget(self.tabs, 1)

        main_tab = QWidget(root)
        main_tab.setProperty("section", "dark")
        main_layout = QVBoxLayout(main_tab)
        main_layout.setContentsMargins(sp(12), sp(10), sp(12), sp(12))
        main_layout.setSpacing(sp(8))

        parameters_tab = QWidget(root)
        parameters_tab.setProperty("section", "dark")
        parameters_layout = QVBoxLayout(parameters_tab)
        parameters_layout.setContentsMargins(sp(12), sp(10), sp(12), sp(12))
        parameters_layout.setSpacing(sp(8))

        results_tab = QWidget(root)
        results_tab.setProperty("section", "dark")
        results_layout = QVBoxLayout(results_tab)
        results_layout.setContentsMargins(sp(12), sp(10), sp(12), sp(12))
        results_layout.setSpacing(sp(8))

        self.tabs.addTab(main_tab, "MAIN")
        self.tabs.addTab(parameters_tab, "PARAMETERS")
        self._results_tab_index = self.tabs.addTab(results_tab, "Residue Table")

        header_row = QHBoxLayout()
        header_row.setContentsMargins(0, 0, 0, 0)
        header_row.setSpacing(sp(8))
        title = QLabel("DaQRotFix", root)
        title.setProperty("role", "title")
        version = QLabel(f"v{_bundle_version()}", root)
        version.setProperty("role", "version-badge")
        header_row.addWidget(title)
        header_row.addStretch(1)
        header_row.addWidget(version, 0, Qt.AlignRight | Qt.AlignVCenter)
        main_layout.addLayout(header_row)

        subtitle = QLabel(
            "Density-assisted Quick Rotamer Fix",
            root,
        )
        subtitle.setProperty("role", "caption")
        main_layout.addWidget(subtitle)

        input_grid = QGridLayout()
        input_grid.setHorizontalSpacing(sp(10))
        input_grid.setVerticalSpacing(sp(8))

        structure_card, structure_layout = self._make_card(
            root, "Structure", sp
        )
        self.structure_menu = QComboBox(root)
        self.structure_menu.setToolTip("Atomic structure containing the rotamers")
        self.structure_menu.currentIndexChanged.connect(self._structure_changed)
        structure_layout.addWidget(self.structure_menu)
        input_grid.addWidget(structure_card, 0, 0)

        map_card, map_layout = self._make_card(root, "Map", sp)
        self.volume_menu = QComboBox(root)
        self.volume_menu.setToolTip("Density map used for side-chain scoring")
        self.volume_menu.currentIndexChanged.connect(self._update_readiness)
        map_layout.addWidget(self.volume_menu)
        input_grid.addWidget(map_card, 0, 1)
        main_layout.addLayout(input_grid)

        mdff_card, mdff_layout = self._make_card(root, "ISOLDE MDFF", sp)
        mdff_row = QHBoxLayout()
        mdff_row.setContentsMargins(0, 0, 0, 0)
        mdff_row.setSpacing(sp(8))
        self.mdff_status = QLabel("Checking model/map association…", root)
        self.mdff_status.setProperty("role", "status-badge")
        self.mdff_status.setProperty("state", "neutral")
        self.mdff_status.setWordWrap(True)
        refresh_button = QPushButton("Refresh", root)
        refresh_button.setProperty("variant", "secondary-gray")
        refresh_button.setToolTip(
            "Recheck the model/map MDFF setup after configuring it in ISOLDE"
        )
        refresh_button.clicked.connect(self._refresh_models)
        mdff_row.addWidget(self.mdff_status, 1)
        mdff_row.addWidget(refresh_button, 0, Qt.AlignRight)
        mdff_layout.addLayout(mdff_row)
        main_layout.addWidget(mdff_card)

        mdff_run_card, mdff_run_layout = self._make_card(
            root, "Local MDFF Optimization", sp
        )
        mdff_tooltip = (
            "<b>Automated Local MDFF</b><br>"
            "Canonical rotamers are tried in decreasing ISOLDE probability. "
            "Each trial receives local MDFF annealing followed by final "
            "minimisation. A candidate "
            "is accepted only if it fixes the target, reduces the local outlier "
            "count and creates no new neighboring outliers.<br><br>"
            "The workflow starts and stops ISOLDE annealing/minimisation runs, "
            "focuses each "
            "active local region, restores rejected coordinates and temporary "
            "backbone restraints, then restores the original view."
        )
        mdff_run_card.setToolTip(mdff_tooltip)

        parameters_header = QLabel("Local MDFF Parameters", root)
        parameters_header.setProperty("role", "title")
        parameters_layout.addWidget(parameters_header)
        parameters_caption = QLabel(
            "These settings are applied to every canonical rotamer trial.", root
        )
        parameters_caption.setProperty("role", "caption")
        parameters_layout.addWidget(parameters_caption)

        parameters_card, parameters_card_layout = self._make_card(
            root, "Optimization Settings", sp
        )
        settings = QGridLayout()
        settings.setHorizontalSpacing(sp(10))
        settings.setVerticalSpacing(sp(10))
        settings.addWidget(QLabel("Local radius (Å)", root), 0, 0)
        self.local_radius = QDoubleSpinBox(root)
        self.local_radius.setRange(2.0, 12.0)
        self.local_radius.setSingleStep(0.5)
        self.local_radius.setValue(DEFAULT_LOCAL_RADIUS)
        settings.addWidget(self.local_radius, 0, 1)

        settings.addWidget(QLabel("Initial temperature (K)", root), 1, 0)
        self.initial_temperature = QDoubleSpinBox(root)
        self.initial_temperature.setDecimals(1)
        self.initial_temperature.setRange(0.0, 300.0)
        self.initial_temperature.setSingleStep(5.0)
        self.initial_temperature.setValue(DEFAULT_INITIAL_TEMPERATURE)
        self.initial_temperature.setToolTip(
            "Temperature used for the first annealing coordinate update."
        )
        settings.addWidget(self.initial_temperature, 1, 1)

        settings.addWidget(QLabel("Annealing updates", root), 2, 0)
        self.annealing_updates = QSpinBox(root)
        self.annealing_updates.setRange(1, 100)
        self.annealing_updates.setSingleStep(5)
        self.annealing_updates.setValue(DEFAULT_ANNEALING_UPDATES)
        self.annealing_updates.setToolTip(
            "Number of ISOLDE dynamics coordinate updates while cooling (1–100)."
        )
        settings.addWidget(self.annealing_updates, 2, 1)

        settings.addWidget(QLabel("Final temperature (K)", root), 3, 0)
        self.final_temperature = QDoubleSpinBox(root)
        self.final_temperature.setDecimals(1)
        self.final_temperature.setRange(0.0, DEFAULT_INITIAL_TEMPERATURE)
        self.final_temperature.setSingleStep(5.0)
        self.final_temperature.setValue(DEFAULT_FINAL_TEMPERATURE)
        self.final_temperature.setToolTip(
            "Temperature reached at the end of the annealing stage."
        )
        self.initial_temperature.valueChanged.connect(
            self.final_temperature.setMaximum
        )
        settings.addWidget(self.final_temperature, 3, 1)

        settings.addWidget(QLabel("Final minimization updates", root), 4, 0)
        self.final_minimization_updates = QSpinBox(root)
        self.final_minimization_updates.setRange(1, 100)
        self.final_minimization_updates.setSingleStep(5)
        self.final_minimization_updates.setValue(
            DEFAULT_FINAL_MINIMIZATION_UPDATES
        )
        self.final_minimization_updates.setToolTip(
            "Number of ISOLDE minimization coordinate updates after cooling (1–100)."
        )
        # Compatibility for external GUI smoke tests written for v0.5.x.
        self.coordinate_updates = self.final_minimization_updates
        settings.addWidget(self.final_minimization_updates, 4, 1)

        settings.addWidget(QLabel("Backbone restraint (kJ mol⁻¹ Å⁻²)", root), 5, 0)
        self.backbone_spring = QDoubleSpinBox(root)
        self.backbone_spring.setRange(10.0, 1000.0)
        self.backbone_spring.setSingleStep(50.0)
        self.backbone_spring.setValue(DEFAULT_BACKBONE_SPRING_CONSTANT / 100.0)
        settings.addWidget(self.backbone_spring, 5, 1)

        settings.addWidget(QLabel("Outlier threshold (%)", root), 6, 0)
        self.outlier_threshold = QDoubleSpinBox(root)
        self.outlier_threshold.setDecimals(3)
        self.outlier_threshold.setRange(0.001, 100.0)
        self.outlier_threshold.setSingleStep(0.01)
        self.outlier_threshold.setValue(DEFAULT_OUTLIER_THRESHOLD * 100.0)
        self.outlier_threshold.setToolTip(
            "A rotamer is an outlier when its ISOLDE P-value is below this percentage"
        )
        self.outlier_threshold.valueChanged.connect(self._refresh_outlier_count)
        settings.addWidget(self.outlier_threshold, 6, 1)
        settings.setColumnStretch(0, 1)
        settings.setColumnStretch(1, 1)
        parameters_card_layout.addLayout(settings)
        parameters_layout.addWidget(parameters_card)

        threshold_note = QLabel(
            "Outliers are residues whose ISOLDE rotamer P-value is below the "
            "selected percentage. The standard ISOLDE cutoff is normally 0.05%.",
            root,
        )
        threshold_note.setProperty("role", "notice")
        threshold_note.setWordWrap(True)
        parameters_layout.addWidget(threshold_note)
        parameters_layout.addStretch(1)

        self.local_mdff_button = QPushButton("Start Optimization", root)
        self.local_mdff_button.setProperty("variant", "primary")
        self.local_mdff_button.setMinimumWidth(sp(160, 130))
        self.local_mdff_button.setToolTip(mdff_tooltip)
        self.local_mdff_button.clicked.connect(self._run_local_mdff)
        mdff_run_layout.addWidget(self.local_mdff_button, 0, Qt.AlignLeft)

        progress_header = QHBoxLayout()
        progress_header.setContentsMargins(0, 0, 0, 0)
        progress_header.setSpacing(sp(8))
        progress_label = QLabel("Progress", root)
        progress_label.setProperty("role", "caption")
        self.outlier_count_label = QLabel("Rotamer outliers: —", root)
        self.outlier_count_label.setProperty("role", "caption")
        progress_header.addWidget(progress_label)
        progress_header.addStretch(1)
        progress_header.addWidget(
            self.outlier_count_label, 0, Qt.AlignRight | Qt.AlignVCenter
        )
        mdff_run_layout.addLayout(progress_header)

        self.progress = QProgressBar(root)
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        mdff_run_layout.addWidget(self.progress)

        self.run_status = QLabel("Ready", root)
        self.run_status.setProperty("role", "caption")
        self.run_status.setWordWrap(True)
        mdff_run_layout.addWidget(self.run_status)
        main_layout.addWidget(mdff_run_card)
        main_layout.addStretch(1)

        results_header = QLabel("Residue Results", root)
        results_header.setProperty("role", "title")
        results_layout.addWidget(results_header)

        self.summary = QLabel("No correction has been run.", root)
        self.summary.setProperty("role", "caption")
        self.summary.setWordWrap(True)
        results_layout.addWidget(self.summary)

        self.results_table = QTableWidget(0, len(self._HEADERS), root)
        self.results_table.setHorizontalHeaderLabels(self._HEADERS)
        self.results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.results_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.results_table.setAlternatingRowColors(True)
        self.results_table.setSortingEnabled(True)
        self.results_table.verticalHeader().setVisible(False)
        self.results_table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        header = self.results_table.horizontalHeader()
        for column in (0, 2):
            header.setSectionResizeMode(column, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.Interactive)
        header.setSectionResizeMode(3, QHeaderView.Interactive)
        header.resizeSection(1, sp(108, 88))
        header.resizeSection(3, sp(130, 108))
        header.setMinimumHeight(sp(48, 40))
        header.setSectionResizeMode(4, QHeaderView.Stretch)
        self.results_table.cellClicked.connect(self._result_clicked)
        results_layout.addWidget(self.results_table, 1)

        detail_card, detail_layout = self._make_card(
            root, "Canonical Target Details", sp
        )
        self.candidate_details = QPlainTextEdit(root)
        self.candidate_details.setReadOnly(True)
        self.candidate_details.setMaximumHeight(sp(170, 110))
        self.candidate_details.setPlaceholderText(
            "Select a residue row to inspect all tested targets."
        )
        detail_layout.addWidget(self.candidate_details)
        results_layout.addWidget(detail_card)

        footer = QHBoxLayout()
        footer.setContentsMargins(sp(2), 0, sp(2), 0)
        hint = QLabel("Select a result row to focus that residue in ChimeraX.", root)
        hint.setProperty("role", "caption")
        footer.addWidget(hint)
        footer.addStretch(1)
        outer.addLayout(footer)

        parent.setLayout(QVBoxLayout())
        parent.layout().setContentsMargins(0, 0, 0, 0)
        parent.layout().addWidget(root)
        self._refresh_models()

    @staticmethod
    def _make_card(root, title, sp):
        card = QFrame(root)
        card.setProperty("card", "dark")
        layout = QVBoxLayout(card)
        layout.setContentsMargins(sp(12), sp(10), sp(12), sp(10))
        layout.setSpacing(sp(6))
        label = QLabel(title, root)
        label.setProperty("role", "field-label")
        layout.addWidget(label)
        return card, layout

    @staticmethod
    def _stylesheet(sp):
        return f"""
            QWidget {{
                background: #000000;
                color: #ffffff;
                font-family: "SF Pro Text", "SF Pro Display", "Helvetica Neue", Helvetica, Arial, sans-serif;
                font-size: 13px;
            }}
            QWidget[section="dark"] {{ background: #000000; color: #ffffff; }}
            QTabWidget::pane {{ border: none; background: transparent; top: 0px; }}
            QTabBar::tab {{
                background: rgba(0, 0, 0, 0.8); color: rgba(255, 255, 255, 0.92);
                border: none; padding: {sp(8, 5)}px {sp(18, 12)}px;
                margin-right: {sp(6, 4)}px; min-width: {sp(92, 76)}px;
                border-radius: 999px; font-size: 12px;
            }}
            QTabBar::tab:selected {{ background: #1d1d1f; color: #ffffff; font-weight: 600; }}
            QTabBar::tab:hover {{ background: rgba(255, 255, 255, 0.10); }}
            QFrame[card="dark"] {{
                background: #1d1d1f; color: #ffffff;
                border: none; border-radius: 10px;
            }}
            QLabel {{ background: transparent; color: #ffffff; }}
            QLabel[role="title"] {{ color: #ffffff; font-size: 20px; font-weight: 600; }}
            QLabel[role="field-label"] {{ color: #ffffff; font-size: 13px; font-weight: 600; }}
            QLabel[role="caption"] {{ font-size: 12px; color: rgba(255,255,255,0.72); }}
            QLabel[role="version-badge"] {{
                color: rgba(255,255,255,0.78); background: rgba(255,255,255,0.10);
                border: 1px solid rgba(255,255,255,0.16); border-radius: 10px;
                padding: 3px 8px; font-size: 12px; font-weight: 600;
            }}
            QLabel[role="status-badge"] {{
                padding: 6px 10px; border-radius: 10px; font-size: 12px; font-weight: 600;
            }}
            QLabel[role="status-badge"][state="ready"] {{
                color: #b6f2c2; background: rgba(52,199,89,0.18);
                border: 1px solid rgba(52,199,89,0.38);
            }}
            QLabel[role="status-badge"][state="warning"] {{
                color: #ffd18a; background: rgba(255,159,10,0.14);
                border: 1px solid rgba(255,159,10,0.38);
            }}
            QLabel[role="status-badge"][state="error"] {{
                color: #ffb4ab; background: rgba(255,69,58,0.14);
                border: 1px solid rgba(255,69,58,0.38);
            }}
            QLabel[role="status-badge"][state="neutral"] {{
                color: rgba(255,255,255,0.72); background: rgba(255,255,255,0.08);
                border: 1px solid rgba(255,255,255,0.14);
            }}
            QLabel[role="notice"] {{
                color: #9ecbff; background: rgba(41,151,255,0.12);
                border: 1px solid rgba(41,151,255,0.30); border-radius: 10px;
                padding: 9px 11px; font-size: 12px;
            }}
            QComboBox, QSpinBox, QDoubleSpinBox {{
                background: #272729; color: #ffffff; border: none;
                border-radius: {sp(12, 8)}px; padding: {sp(6, 4)}px {sp(10, 7)}px;
                min-height: {sp(22, 18)}px; selection-background-color: #0071e3;
            }}
            QSpinBox QLineEdit, QDoubleSpinBox QLineEdit {{
                background: transparent; color: #ffffff;
                selection-background-color: #0071e3; selection-color: #ffffff;
            }}
            QComboBox QAbstractItemView {{
                background: #272729; color: #ffffff; border: 1px solid #4a4a4d;
                selection-background-color: #0071e3; selection-color: #ffffff;
            }}
            QComboBox::drop-down {{ border: 0px; width: 24px; }}
            QPushButton {{
                background: #ffffff; color: #1d1d1f;
                border: 1px solid rgba(0,0,0,0.08); border-radius: {sp(14, 9)}px;
                padding: {sp(6, 4)}px {sp(13, 8)}px; min-height: {sp(22, 18)}px;
            }}
            QPushButton:hover {{ background: #fafafc; }}
            QPushButton[variant="primary"] {{
                background: #0071e3; color: #ffffff; border: 1px solid transparent;
            }}
            QPushButton[variant="primary"]:hover {{ background: #0077ed; }}
            QPushButton[variant="secondary-gray"] {{
                background: #3a3a3c; color: #ffffff;
                border: 1px solid rgba(255,255,255,0.14);
            }}
            QPushButton[variant="secondary-gray"]:hover {{ background: #4a4a4d; }}
            QPushButton:disabled {{
                background: #2a2a2d; color: rgba(255,255,255,0.42);
                border: 1px solid rgba(255,255,255,0.08);
            }}
            QProgressBar {{
                background: #111113; color: #ffffff; border: 1px solid #2f2f32;
                border-radius: 8px; min-height: 16px; text-align: center; font-size: 11px;
            }}
            QProgressBar::chunk {{ background: #0071e3; border-radius: 7px; }}
            QTableWidget {{
                background: #111113; color: #ffffff; border: 1px solid #2f2f32;
                border-radius: 12px; gridline-color: #3a3a3d;
                alternate-background-color: #1a1a1c;
                selection-background-color: #2f6fed; selection-color: #ffffff;
            }}
            QTableWidget::item {{ padding: 6px 8px; border: none; }}
            QHeaderView::section {{
                background: #1d1d1f; color: #ffffff; border: none;
                border-right: 1px solid #3a3a3d; border-bottom: 1px solid #3a3a3d;
                padding: 8px 10px; font-size: 13px; font-weight: 600;
            }}
            QTableCornerButton::section {{
                background: #1d1d1f; border: none;
                border-right: 1px solid #3a3a3d; border-bottom: 1px solid #3a3a3d;
            }}
            QPlainTextEdit {{
                background: #111113; color: #ffffff; border: 1px solid #2f2f32;
                border-radius: 10px; padding: 7px; selection-background-color: #2f6fed;
                font-family: Menlo, Monaco, monospace; font-size: 12px;
            }}
            QToolTip {{ background: #272729; color: #ffffff; border: 1px solid #4a4a4d; }}
        """

    def _on_models_changed(self, trigger_name, data):
        self._refresh_models()

    def _models(self):
        from chimerax.atomic import AtomicStructure
        from chimerax.map import Volume

        models = self.session.models.list()
        return (
            [model for model in models if isinstance(model, AtomicStructure)],
            [model for model in models if isinstance(model, Volume)],
        )

    @staticmethod
    def _model_label(model):
        return f"#{model.id_string} {model.name}"

    def _refresh_models(self):
        current_structure = self.structure_menu.currentData()
        current_volume = self.volume_menu.currentData()
        structures, volumes = self._models()

        self.structure_menu.blockSignals(True)
        self.volume_menu.blockSignals(True)
        self.structure_menu.clear()
        self.volume_menu.clear()
        for structure in structures:
            self.structure_menu.addItem(self._model_label(structure), structure)
        for volume in volumes:
            self.volume_menu.addItem(self._model_label(volume), volume)
        self._restore_choice(self.structure_menu, current_structure)
        self._restore_choice(self.volume_menu, current_volume)
        self.structure_menu.blockSignals(False)
        self.volume_menu.blockSignals(False)
        self._choose_matching_map()
        self._initialize_outlier_threshold()
        self._update_readiness()
        self._refresh_outlier_count()

    def _initialize_outlier_threshold(self):
        if self._threshold_initialized:
            return
        try:
            from chimerax.isolde import session_extensions

            rotamer_mgr = session_extensions.get_rotamer_mgr(self.session)
            cutoff = float(rotamer_mgr.cutoffs[1])
        except Exception:
            return
        if 0.0 < cutoff <= 1.0:
            self.outlier_threshold.setValue(cutoff * 100.0)
            self._threshold_initialized = True

    @staticmethod
    def _restore_choice(menu, model):
        if model is None:
            return
        for index in range(menu.count()):
            if menu.itemData(index) is model:
                menu.setCurrentIndex(index)
                return

    def _structure_changed(self):
        self._choose_matching_map()
        self._update_readiness()
        self._refresh_outlier_count()

    def _set_outlier_count(self, count):
        value = "—" if count is None else str(int(count))
        self.outlier_count_label.setText(f"Rotamer outliers: {value}")

    def _refresh_outlier_count(self, *_args, force=False):
        if self._running and not force:
            return
        structure = self.structure_menu.currentData()
        if structure is None:
            self._set_outlier_count(None)
            return
        try:
            from chimerax.isolde import session_extensions

            rotamer_mgr = session_extensions.get_rotamer_mgr(self.session)
            threshold = self.outlier_threshold.value() / 100.0
            count = len(
                _outliers_in_residues(
                    rotamer_mgr, structure.residues, threshold
                )
            )
        except Exception:
            self._set_outlier_count(None)
            return
        self._set_outlier_count(count)

    def _choose_matching_map(self):
        structure = self.structure_menu.currentData()
        if structure is None:
            return
        from chimerax.isolde import session_extensions

        current = self.volume_menu.currentData()
        if current is not None and session_extensions.get_mdff_mgr(
            structure, current, create=False
        ) is not None:
            return
        for index in range(self.volume_menu.count()):
            volume = self.volume_menu.itemData(index)
            if session_extensions.get_mdff_mgr(
                structure, volume, create=False
            ) is not None:
                self.volume_menu.setCurrentIndex(index)
                return

    def _current_mdff_mgr(self):
        structure = self.structure_menu.currentData()
        volume = self.volume_menu.currentData()
        if structure is None or volume is None:
            return None
        from chimerax.isolde import session_extensions

        return session_extensions.get_mdff_mgr(structure, volume, create=False)

    def _set_status_badge(self, text, state):
        self.mdff_status.setText(text)
        self.mdff_status.setProperty("state", state)
        self.mdff_status.style().unpolish(self.mdff_status)
        self.mdff_status.style().polish(self.mdff_status)

    def _update_readiness(self):
        structure = self.structure_menu.currentData()
        volume = self.volume_menu.currentData()
        isolde = getattr(self.session, "isolde", None)

        if structure is None or volume is None:
            self._set_status_badge(
                "Open ISOLDE, select a structure and density map, and set up MDFF.",
                "warning",
            )
            local_ready = False
        elif isolde is not None and isolde.simulation_running:
            self._set_status_badge(
                "An ISOLDE simulation is active. Stop it before changing coordinates.",
                "error",
            )
            local_ready = False
        elif self._current_mdff_mgr() is None:
            self._set_status_badge(
                "No MDFF setup for this model/map. Set it up in ISOLDE, then "
                "click Refresh.",
                "error",
            )
            local_ready = False
        elif isolde is None:
            self._set_status_badge(
                "Open ISOLDE before starting automated local MDFF.", "warning"
            )
            local_ready = False
        elif not bool(getattr(self._current_mdff_mgr(), "enabled", False)):
            self._set_status_badge(
                "Enable this model/map's MDFF potential in ISOLDE, then click "
                "Refresh.",
                "warning",
            )
            local_ready = False
        else:
            self._set_status_badge(
                "ISOLDE MDFF setup found and enabled. Ready for local optimization.",
                "ready",
            )
            local_ready = True
        self.local_mdff_button.setEnabled(local_ready and not self._running)

    def _local_mdff_progress(self, stage, completed, total, residue, detail, extra):
        self.progress.setRange(0, max(1, total))
        self.progress.setValue(completed)
        if stage == "starting":
            self._set_outlier_count(total)
            self.run_status.setText(
                f"Found {total} original outlier(s). Preparing local MDFF…"
            )
        elif stage == "candidate":
            self.run_status.setText(
                f"{residue_label(residue)}: local MDFF candidate "
                f"{detail}/{extra} ({completed + 1}/{total})"
            )
        elif stage == "residue":
            self._refresh_outlier_count(force=True)
            status = "accepted" if detail.corrected else "unchanged"
            self.run_status.setText(
                f"{residue_label(residue)}: {status} ({completed}/{total})"
            )
        elif stage == "complete":
            self._refresh_outlier_count(force=True)
            self.run_status.setText("Automated local MDFF optimization complete.")
        QApplication.processEvents()

    def _run_local_mdff(self):
        structure = self.structure_menu.currentData()
        volume = self.volume_menu.currentData()
        if structure is None or volume is None:
            QMessageBox.warning(
                self.tool_window.ui_area,
                "DaQRotFix",
                "Choose both a structure and a density map.",
            )
            return

        answer = QMessageBox.question(
            self.tool_window.ui_area,
            "Automated Local MDFF",
            "This action will automatically start and stop short ISOLDE "
            "annealing and minimisation runs for each rotamer candidate. Continue?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        )
        if answer != QMessageBox.Yes:
            return

        self._running = True
        initial_view = _capture_view(self.session)
        self.local_mdff_button.setEnabled(False)
        self.results_table.setRowCount(0)
        self.candidate_details.clear()
        self.summary.setText("Automated local MDFF optimization in progress…")
        self.progress.setRange(0, 0)
        self.outlier_count_label.setText("Rotamer outliers: scanning…")
        self.run_status.setText("Finding original rotamer outliers…")
        QApplication.processEvents()
        try:
            try:
                report = optimize_sidechains_with_local_mdff(
                    self.session,
                    structure,
                    volume,
                    local_radius=self.local_radius.value(),
                    initial_temperature=self.initial_temperature.value(),
                    annealing_updates=self.annealing_updates.value(),
                    final_temperature=self.final_temperature.value(),
                    final_minimization_updates=(
                        self.final_minimization_updates.value()
                    ),
                    # GUI is expressed per A^2; ISOLDE stores per nm^2.
                    backbone_spring_constant=self.backbone_spring.value() * 100.0,
                    # GUI is percentage; the validator returns a 0..1 P-value.
                    outlier_threshold=self.outlier_threshold.value() / 100.0,
                    progress_callback=self._local_mdff_progress,
                )
            finally:
                _restore_view(self.session, initial_view)
            self._report = report
            self._populate_local_mdff_results(report)
            self._log_local_mdff_report(report)
            self.tabs.setCurrentIndex(self._results_tab_index)
        except (UserError, RuntimeError, ValueError) as error:
            self.run_status.setText("Automated local MDFF optimization failed.")
            self.summary.setText(str(error))
            self.session.logger.error(str(error))
            QMessageBox.critical(
                self.tool_window.ui_area, "DaQRotFix", str(error)
            )
        finally:
            self._running = False
            self._update_readiness()
            self._refresh_outlier_count()

    def _populate_local_mdff_results(self, report):
        self.results_table.setHorizontalHeaderLabels(
            self._HEADERS
        )
        self.results_table.setSortingEnabled(False)
        self.results_table.setRowCount(len(report.residues))
        for row, result in enumerate(report.residues):
            selected = result.selected
            if result.resolved_before_trial:
                values = (
                    residue_label(result.residue),
                    "—",
                    "—",
                    str(result.baseline_local_outlier_count),
                    "Already resolved by an earlier local optimization",
                )
            elif selected is None:
                values = (
                    residue_label(result.residue),
                    "—",
                    "—",
                    f"{result.baseline_local_outlier_count} → unchanged",
                    "Unchanged — no safe local improvement",
                )
            else:
                values = (
                    residue_label(result.residue),
                    selected.target_name,
                    f"{selected.probability:.3g}%",
                    f"{result.baseline_local_outlier_count} → "
                    f"{selected.local_outlier_count}",
                    "Accepted — target fixed, no new neighboring outliers",
                )
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setData(Qt.UserRole, result)
                if column in (2, 3):
                    item.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
                self.results_table.setItem(row, column, item)
        self.results_table.setSortingEnabled(True)
        self.summary.setText(
            f"Optimized {report.corrected_count} of "
            f"{report.original_outlier_count} original rotamer outlier(s) "
            f"using probability-ordered local ISOLDE MDFF trials "
            f"(outlier threshold {report.outlier_threshold * 100.0:.3g}%)."
        )
        if report.residues:
            self.results_table.selectRow(0)
            self._show_local_mdff_details(report.residues[0])
        else:
            self.candidate_details.setPlainText("No original rotamer outliers were found.")

    def _log_local_mdff_report(self, report):
        logger = self.session.logger
        logger.info(
            f"DaQRotFix: found {report.original_outlier_count} "
            f"original rotamer outlier(s) at threshold "
            f"{report.outlier_threshold * 100.0:.3g}%."
        )
        for result in report.residues:
            selected = result.selected
            if result.resolved_before_trial:
                logger.info(
                    f"{residue_label(result.residue)}: already resolved by an "
                    "earlier local optimization"
                )
            elif selected is None:
                logger.warning(
                    f"{residue_label(result.residue)}: unchanged; no candidate "
                    "fixed the target without creating a neighboring outlier"
                )
            else:
                logger.info(
                    f"{residue_label(result.residue)}: selected "
                    f"{selected.target_name} ({selected.probability:.3g}%), "
                    f"local outliers {result.baseline_local_outlier_count} -> "
                    f"{selected.local_outlier_count}"
                )

    def _result_clicked(self, row, column):
        item = self.results_table.item(row, 0)
        if item is None:
            return
        result = item.data(Qt.UserRole)
        if result is None:
            return
        self._show_local_mdff_details(result)
        residue = result.residue
        if getattr(residue, "deleted", False):
            return
        self.session.selection.clear()
        residue.atoms.selected = True
        try:
            from chimerax.core.commands import run

            run(self.session, "view sel", log=False)
        except UserError:
            pass

    def _show_local_mdff_details(self, result):
        if result.resolved_before_trial:
            self.candidate_details.setPlainText(
                "Already resolved by an earlier local optimization; no candidates tested."
            )
            return
        selected_index = (
            result.selected.target_index if result.selected is not None else None
        )
        lines = []
        for candidate in result.candidates:
            marker = "★" if candidate.target_index == selected_index else " "
            if candidate.completed:
                outliers = ", ".join(candidate.outlier_labels) or "none"
                decision = (
                    "accepted"
                    if candidate.selection_reason is None
                    else f"rejected — {candidate.selection_reason}"
                )
                lines.append(
                    f"{marker} {candidate.target_name} "
                    f"({candidate.probability:.3g}%): "
                    f"local outliers {result.baseline_local_outlier_count} → "
                    f"{candidate.local_outlier_count} [{outliers}]; {decision}"
                )
            else:
                lines.append(
                    f"{marker} {candidate.target_name} "
                    f"({candidate.probability:.3g}%): failed — "
                    f"{candidate.rejection_reason}"
                )
        self.candidate_details.setPlainText("\n".join(lines))
