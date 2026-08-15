"""Manual GUI smoke test for parameter-label contrast."""

from Qt.QtGui import QPalette
from Qt.QtWidgets import QApplication, QLabel
from chimerax.core.commands import run


run(session, "ui tool show DaQRotFix")
QApplication.processEvents()
tool = session.tools.find_by_class(
    __import__(
        "chimerax.daq_rot_fix.tool", fromlist=["DaQRotFixTool"]
    ).DaQRotFixTool
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
assert tool.final_minimization_updates.value() == 5
assert tool._HEADERS[1] == "Selected\nRotamer"
assert tool._HEADERS[3] == "Local Outliers\n(Before → After)"
assert tool.results_table.columnWidth(1) < 130
assert tool.results_table.columnWidth(3) < 150
assert "ISOLDE" in tool.mdff_status.text()
session.logger.info(f"DAQROTFIX_PARAMETER_LABELS_WHITE {found}")
