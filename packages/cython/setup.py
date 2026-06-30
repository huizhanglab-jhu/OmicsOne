# setup.py
from pathlib import Path
from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from Cython.Build import cythonize
import numpy as np
import os

# ── 基本路径（仅 THIS_DIR 用绝对路径；所有传给 setuptools/cython 的路径统一转相对 + 正斜杠） ──
THIS_DIR = Path(__file__).parent.resolve()         # .../packages/cython
SRC_DIR  = (THIS_DIR / "../../src").resolve()      # 项目源码根（发现文件时用）

def to_rel(p) -> str:
    """把任何路径转成相对 setup.py 的路径，且用 '/' 分隔。"""
    return os.path.relpath(str(p), start=THIS_DIR).replace("\\", "/")

# ── 需要“像 utils 一样处理”的包名列表：这些包内除 __init__.py 外的 .py 都会被 cythonize，并从纯 .py 构建中剔除 ──
TARGET_PKG_NAMES = [
    "omicsone_streamlit.utils",
    "omicsone_streamlit.myviews",   # ← 新增：像 utils 一样处理 myviews
    "omicsone_streamlit.mypages",  # ← 新增：像 utils 一样处理 mypages
    "omicsone_streamlit.plots",
    # 以后再加别的：比如 "omicsone_streamlit.plots", "omicsone_streamlit.services", ...
]

# 把包名转成具体目录路径
def pkg_to_dir(pkg_name: str) -> Path:
    return SRC_DIR / Path(pkg_name.replace(".", "/"))

# ── 收集要 cythonize 的源文件（各目标包里，除 __init__.py 外的所有 .py） ──
to_cythonize = []
for pkg in TARGET_PKG_NAMES:
    pkg_dir = pkg_to_dir(pkg)
    for p in pkg_dir.rglob("*.py"):
        if p.name == "__init__.py":
            continue
        to_cythonize.append(to_rel(p))

if not to_cythonize:
    print("[setup.py] No .py sources to cythonize in target packages (excluding __init__.py).")

# ── 运行 cythonize：使用相对 build_dir，避免把绝对路径写进 Extension ──
ext_modules = cythonize(
    to_cythonize,
    compiler_directives={"language_level": "3", "embedsignature": True},
    build_dir=to_rel(THIS_DIR / "build" / "cython"),  # 生成的 .c 放到构建目录
    force=True,                                       # 避免复用历史产物
)

# 保险：把 cythonize 生成的 Extension 中的 sources/depends 统一为相对 + 正斜杠
for ext in ext_modules:
    if getattr(ext, "sources", None):
        ext.sources = [to_rel(s) for s in ext.sources]
    if getattr(ext, "depends", None):
        ext.depends = [to_rel(d) for d in ext.depends]

# ── 覆盖 build_py：在构建纯 Python 模块阶段，剔除目标包中除 __init__.py 外的 .py，确保 wheel 里只留 .pyd ──
class build_py_exclude_targets(_build_py):
    """在 build 阶段剔除 TARGET_PKG_NAMES 中（除 __init__.py 外）的纯 .py 源。"""
    def find_package_modules(self, package, package_dir):
        modules = super().find_package_modules(package, package_dir)
        filtered = []
        for pkg, mod, fname in modules:
            # 例如：pkg == "omicsone_streamlit.utils" 或 "omicsone_streamlit.myviews"
            if pkg in TARGET_PKG_NAMES and mod != "__init__":
                # 跳过该 .py，让同名 .pyd 生效
                continue
            filtered.append((pkg, mod, fname))
        return filtered

# ── setup ──
setup(
    # 其余元数据（name/version 等）建议放在 pyproject.toml
    ext_modules=ext_modules,
    include_dirs=[np.get_include()],
    cmdclass={"build_py": build_py_exclude_targets},
    # zip_safe=False,  # 如有需要可开启
)
