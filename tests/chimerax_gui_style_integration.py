"""Manual GUI smoke test for parameter-label contrast."""

from Qt.QtGui import QPalette
from Qt.QtWidgets import QApplication, QLabel
from chimerax.core.commands import run


run(session, "ui tool show Rotamer Optimizer")
QApplication.processEvents()
tool = session.tools.find_by_class(
    __import__(
        "chimerax.correct_rotamer.tool", fromlist=["CorrectRotamerTool"]
    ).CorrectRotamerTool
)[0]

expected = {
    "Local radius (Å)",
    "Initial temperature (K)",
    "Annealing updates",
    "Final temperature (K)",
    "Final minimization updates",
    "Backbone restraint (kJ mol⁻¹ Å⁻²)",
    "Outlier threshold (%)",
}
found = {}
for label in tool.tool_window.ui_area.findChildren(QLabel):
    if label.text() in expected:
        found[label.text()] = label.palette().color(QPalette.WindowText).name()

assert set(found) == expected, found
assert all(color.lower() == "#ffffff" for color in found.values()), found
visible_text = {
    label.text()
    for label in tool.tool_window.ui_area.findChildren(QLabel)
    if label.isVisible()
}
assert not any(text.startswith("Tries canonical rotamers") for text in visible_text)
assert not any(text.startswith("This workflow explicitly") for text in visible_text)
assert "Canonical rotamers are tried" in tool.local_mdff_button.toolTip()
assert tool.initial_temperature.value() == 20.0
assert tool.annealing_updates.minimum() == 1
assert tool.annealing_updates.maximum() == 100
assert tool.annealing_updates.value() == 20
assert tool.final_temperature.value() == 5.0
assert tool.final_temperature.maximum() == 20.0
assert tool.final_minimization_updates.minimum() == 1
assert tool.final_minimization_updates.maximum() == 100
assert tool.final_minimization_updates.value() == 20
session.logger.info(f"CORRECTROTAMER_PARAMETER_LABELS_WHITE {found}")
