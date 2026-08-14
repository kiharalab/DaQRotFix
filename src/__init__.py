"""ChimeraX CorrectRotamer bundle."""

from chimerax.core.toolshed import BundleAPI


class _CorrectRotamerAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def register_command(bundle_info, command_info, logger):
        from chimerax.core.commands import register
        from .cmd import (
            correct_rotamers_cmd,
            correct_rotamers_desc,
            correct_rotamers_mdff_cmd,
            correct_rotamers_mdff_desc,
        )

        commands = {
            "correctrotamers": (correct_rotamers_desc, correct_rotamers_cmd),
            "correctrotamers mdff": (
                correct_rotamers_mdff_desc,
                correct_rotamers_mdff_cmd,
            ),
        }
        try:
            description, function = commands[command_info.name]
        except KeyError:
            raise ValueError(f"Unknown command: {command_info.name}")
        if description.synopsis is None:
            description.synopsis = command_info.synopsis
        register(
            command_info.name,
            description,
            function,
            logger=logger,
        )

    @staticmethod
    def start_tool(session, bundle_info, tool_info):
        if tool_info.name not in ("Rotamer Optimizer", "Correct Rotamers"):
            raise ValueError(f"Unknown tool: {tool_info.name}")
        from .tool import CorrectRotamerTool

        return CorrectRotamerTool(session, tool_info.name)


bundle_api = _CorrectRotamerAPI()
