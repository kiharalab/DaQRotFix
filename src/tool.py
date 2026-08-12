"""DAQplugin-styled ChimeraX graphical front end for CorrectRotamer."""

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
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .cmd import _log_report
from .correction import correct_rotamer_outliers, residue_label


def _bundle_version():
    if package_version is None:
        return "0.2.0"
    try:
        return package_version("ChimeraX-CorrectRotamer")
    except Exception:
        return "0.2.0"


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


class CorrectRotamerTool(ToolInstance):
    SESSION_ENDURING = False
    SESSION_SAVE = False
    HELP = "help:command:correctrotamers"

    _HEADERS = (
        "Residue",
        "Selected Rotamer",
        "Mean Density",
        "Valid / Tested",
        "Result",
    )

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "CorrectRotamer"
        self._running = False
        self._report = None

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

        results_tab = QWidget(root)
        results_tab.setProperty("section", "dark")
        results_layout = QVBoxLayout(results_tab)
        results_layout.setContentsMargins(sp(12), sp(10), sp(12), sp(12))
        results_layout.setSpacing(sp(8))

        self.tabs.addTab(main_tab, "MAIN")
        self._results_tab_index = self.tabs.addTab(results_tab, "Residue Table")

        header_row = QHBoxLayout()
        header_row.setContentsMargins(0, 0, 0, 0)
        header_row.setSpacing(sp(8))
        title = QLabel("CorrectRotamer", root)
        title.setProperty("role", "title")
        version = QLabel(f"v{_bundle_version()}", root)
        version.setProperty("role", "version-badge")
        header_row.addWidget(title)
        header_row.addStretch(1)
        header_row.addWidget(version, 0, Qt.AlignRight | Qt.AlignVCenter)
        main_layout.addLayout(header_row)

        subtitle = QLabel(
            "Density-guided canonical rotamer selection using ISOLDE validation",
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
        refresh_button.setToolTip("Refresh open structures and maps")
        refresh_button.clicked.connect(self._refresh_models)
        mdff_row.addWidget(self.mdff_status, 1)
        mdff_row.addWidget(refresh_button, 0, Qt.AlignRight)
        mdff_layout.addLayout(mdff_row)
        main_layout.addWidget(mdff_card)

        run_card, run_layout = self._make_card(root, "Rotamer Correction", sp)
        explanation = QLabel(
            "Tests every canonical target from the original coordinates and "
            "selects the valid target with the highest mean side-chain density.",
            root,
        )
        explanation.setWordWrap(True)
        explanation.setProperty("role", "caption")
        run_layout.addWidget(explanation)

        self.run_button = QPushButton("Correct Rotamer Outliers", root)
        self.run_button.setProperty("variant", "primary")
        self.run_button.setMinimumWidth(sp(220, 160))
        self.run_button.clicked.connect(self._run)
        run_layout.addWidget(self.run_button, 0, Qt.AlignLeft)

        self.progress = QProgressBar(root)
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        run_layout.addWidget(self.progress)

        self.run_status = QLabel("Ready", root)
        self.run_status.setProperty("role", "caption")
        self.run_status.setWordWrap(True)
        run_layout.addWidget(self.run_status)
        main_layout.addWidget(run_card)

        notice = QLabel(
            "CorrectRotamer does not start, pause, or minimize an ISOLDE "
            "simulation. Final refinement remains under your control in ISOLDE.",
            root,
        )
        notice.setProperty("role", "notice")
        notice.setWordWrap(True)
        main_layout.addWidget(notice)
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
        for column in range(4):
            header.setSectionResizeMode(column, QHeaderView.ResizeToContents)
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
            QFrame[card="dark"] {{ background: #1d1d1f; border: none; border-radius: 10px; }}
            QLabel {{ background: transparent; color: inherit; }}
            QLabel[role="title"] {{ font-size: 20px; font-weight: 600; }}
            QLabel[role="field-label"] {{ font-size: 13px; font-weight: 600; }}
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
            QComboBox {{
                background: #272729; color: #ffffff; border: none;
                border-radius: {sp(12, 8)}px; padding: {sp(6, 4)}px {sp(10, 7)}px;
                min-height: {sp(22, 18)}px; selection-background-color: #0071e3;
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
        self._update_readiness()

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
                "Open and select both a structure and a density map.", "warning"
            )
            ready = False
        elif isolde is not None and isolde.simulation_running:
            self._set_status_badge(
                "An ISOLDE simulation is active. Stop it before changing coordinates.",
                "error",
            )
            ready = False
        elif self._current_mdff_mgr() is None:
            self._set_status_badge(
                "No existing MDFF manager for this model/map pair.", "error"
            )
            ready = False
        else:
            self._set_status_badge(
                "Existing MDFF manager found. Ready for density scoring.", "ready"
            )
            ready = True
        self.run_button.setEnabled(ready and not self._running)

    def _progress_update(self, stage, completed, total, residue, result):
        self.progress.setRange(0, max(1, total))
        self.progress.setValue(completed)
        if stage == "enumerating":
            if residue is None:
                self.run_status.setText(
                    f"Found {total} original outlier(s). Testing canonical targets…"
                )
            else:
                action = "corrected" if result.corrected else "unchanged"
                self.run_status.setText(
                    f"{residue_label(residue)}: {action} ({completed}/{total})"
                )
        elif stage == "complete":
            self.run_status.setText(
                "Correction complete. ISOLDE simulation was not started."
            )
        QApplication.processEvents()

    def _run(self):
        structure = self.structure_menu.currentData()
        volume = self.volume_menu.currentData()
        if structure is None or volume is None:
            QMessageBox.warning(
                self.tool_window.ui_area,
                "CorrectRotamer",
                "Choose both a structure and a density map.",
            )
            return

        self._running = True
        self.run_button.setEnabled(False)
        self.results_table.setRowCount(0)
        self.candidate_details.clear()
        self.summary.setText("Correction in progress…")
        self.progress.setRange(0, 0)
        self.run_status.setText("Finding original rotamer outliers…")
        QApplication.processEvents()
        try:
            report = correct_rotamer_outliers(
                self.session,
                structure,
                volume,
                progress_callback=self._progress_update,
            )
            self._report = report
            self._populate_results(report)
            _log_report(self.session, report)
            self.tabs.setCurrentIndex(self._results_tab_index)
        except (UserError, RuntimeError, ValueError) as error:
            self.run_status.setText("Correction failed.")
            self.summary.setText(str(error))
            self.session.logger.error(str(error))
            QMessageBox.critical(
                self.tool_window.ui_area, "CorrectRotamer", str(error)
            )
        finally:
            self._running = False
            self._update_readiness()

    def _populate_results(self, report):
        self.results_table.setSortingEnabled(False)
        self.results_table.setRowCount(len(report.residues))
        for row, result in enumerate(report.residues):
            selected = result.selected
            valid_count = sum(candidate.valid for candidate in result.candidates)
            values = (
                residue_label(result.residue),
                selected.target_name if selected else "—",
                f"{selected.density_score:.6g}" if selected else "—",
                f"{valid_count} / {len(result.candidates)}",
                "Corrected" if selected else "Unchanged — no valid in-map target",
            )
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setData(Qt.UserRole, result)
                if column in (2, 3):
                    item.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
                self.results_table.setItem(row, column, item)
        self.results_table.setSortingEnabled(True)
        self.summary.setText(
            f"Corrected {report.corrected_count} of "
            f"{report.original_outlier_count} original rotamer outlier(s). "
            "ISOLDE was not started or modified."
        )
        if report.residues:
            self.results_table.selectRow(0)
            self._show_candidate_details(report.residues[0])
        else:
            self.candidate_details.setPlainText("No original rotamer outliers were found.")

    def _result_clicked(self, row, column):
        item = self.results_table.item(row, 0)
        if item is None:
            return
        result = item.data(Qt.UserRole)
        if result is None:
            return
        self._show_candidate_details(result)
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

    def _show_candidate_details(self, result):
        selected_index = (
            result.selected.target_index if result.selected is not None else None
        )
        lines = []
        for candidate in result.candidates:
            marker = "★" if candidate.target_index == selected_index else " "
            if candidate.valid:
                lines.append(
                    f"{marker} {candidate.target_name}: "
                    f"mean density {candidate.density_score:.6g}"
                )
            else:
                lines.append(
                    f"{marker} {candidate.target_name}: rejected — "
                    f"{candidate.rejection_reason}"
                )
        self.candidate_details.setPlainText("\n".join(lines))
