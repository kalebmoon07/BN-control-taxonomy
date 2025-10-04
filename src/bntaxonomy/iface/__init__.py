
import glob
import os
from bntaxonomy.utils.log import main_logger

def load_tools():
    basedirs = [os.path.dirname(__file__)]
    __TOOLS.clear()
    for basedir in basedirs:
        for toolfile in glob.glob(f"{basedir}/*.py"):
            toolname = os.path.basename(toolfile)[:-3]
            if toolname[0] == "_":
                continue
            try:
                __import__(f"bntaxonomy.iface.{toolname}")
            except Exception as e:
                main_logger.warning(f"Fail to load tool '{toolname}' ({e})")

__TOOLS = []

def tool_names():
    return [t.name for t in __TOOLS]

def registered_tools():
    return iter(__TOOLS)

def register_tool(toolcls):
    if not hasattr(toolcls, "uses_cache"):
        toolcls.uses_cache = False
    if not hasattr(toolcls, "name"):
        toolcls.name = toolcls.__name__
    __TOOLS.append(toolcls)
    main_logger.info(f"Registered tool {toolcls.name}")
    return toolcls
