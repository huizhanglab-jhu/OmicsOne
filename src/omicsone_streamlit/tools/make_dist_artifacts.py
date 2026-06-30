# tools/make_dist_artifacts.py
import shutil, os, glob, pathlib

SRC = pathlib.Path(".")
DST = pathlib.Path("dist_artifacts")
DST.mkdir(exist_ok=True)

# 总是完整拷贝的目录/文件（UI层）
ALWAYS_COPY = ["app.py", "config", "mypages", "settings.toml", "README.md", "LICENSE"]

# 只拷贝编译产物的目录
CORE_DIRS = ["utils", "services", "components", "myviews", "plots"]

def copy_ui():
    for p in ALWAYS_COPY:
        src = SRC / p
        dst = DST / p
        if src.is_dir():
            shutil.copytree(src, dst, dirs_exist_ok=True)
        elif src.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)

def copy_core():
    for d in CORE_DIRS:
        src_dir = SRC / d
        if not src_dir.exists(): continue
        for root, _, files in os.walk(src_dir):
            rootp = pathlib.Path(root)
            rel = rootp.relative_to(SRC)
            (DST / rel).mkdir(parents=True, exist_ok=True)
            # 保留 __init__.py
            if (rootp / "__init__.py").exists():
                shutil.copy2(rootp / "__init__.py", DST / rel / "__init__.py")
            # 拷贝 .pyd/.so
            for ext in ("*.pyd", "*.so", "*.dll"):
                for f in glob.glob(str(rootp / ext)):
                    shutil.copy2(f, DST / rel / pathlib.Path(f).name)

if __name__ == "__main__":
    if DST.exists(): shutil.rmtree(DST)
    DST.mkdir()
    copy_ui()
    copy_core()
    print("OK -> dist_artifacts/")
