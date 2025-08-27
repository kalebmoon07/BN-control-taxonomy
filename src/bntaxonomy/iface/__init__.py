
import glob
import os
from bntaxonomy.utils.log import main_logger

def load_tools(extra_dirs=[]):
    basedirs = [os.path.dirname(__file__)] + extra_dirs
    for basedir in basedirs:
        for toolfile in glob.glob(f"{basedir}/*.py"):
            if toolfile[0] == "_":
                continue
            try:
                __import__(f".{toolfile[:-3]}")
            except Exception as e:
                main_logger.warning(f"Fail to load tool file {toolfile} ({e})")

registered_tools = []

def register_tool(toolcls):
    registered_tools.append(toolcls)
    return toolcls
