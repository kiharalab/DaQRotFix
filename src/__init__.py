"""ChimeraX CorrectRotamer bundle."""

from chimerax.core.toolshed import BundleAPI


class _CorrectRotamerAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def register_command(bundle_info, command_info, logger):
        from chimerax.core.commands import register
        from .cmd import correct_rotamers_cmd, correct_rotamers_desc

        if command_info.name != "correctrotamers":
            raise ValueError(f"Unknown command: {command_info.name}")
        if correct_rotamers_desc.synopsis is None:
            correct_rotamers_desc.synopsis = command_info.synopsis
        register(
            command_info.name,
            correct_rotamers_desc,
            correct_rotamers_cmd,
            logger=logger,
        )

    @staticmethod
    def start_tool(session, bundle_info, tool_info):
        if tool_info.name != "Correct Rotamers":
            raise ValueError(f"Unknown tool: {tool_info.name}")
        from .tool import CorrectRotamerTool

        return CorrectRotamerTool(session, tool_info.name)


bundle_api = _CorrectRotamerAPI()
