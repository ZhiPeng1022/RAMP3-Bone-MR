# 解压 eBMD Morris 2018 tar.gz 并列出内部文件

import gzip
import os
import tarfile


TAR = os.path.join(r"data/derived", "Morrisetal2018.NatGen.SumStats.tar.gz")
OUT_DIR = os.path.join(r"data/derived", "ebmd_unpacked")


def main():
    if not os.path.exists(TAR) or os.path.getsize(TAR) < 1062749253:
        print("tar 文件未完整，当前大小：", os.path.getsize(TAR) if os.path.exists(TAR) else "不存在")
        return

    os.makedirs(OUT_DIR, exist_ok=True)
    with tarfile.open(TAR, "r:gz") as tf:
        members = tf.getmembers()
        print("tar 成员数：", len(members))
        for m in members[:20]:
            print(m.name, m.size)
        tf.extractall(OUT_DIR)

    print("已解压到：", OUT_DIR)
    # 查看解压后的文件
    for root, _, files in os.walk(OUT_DIR):
        for f in files[:20]:
            path = os.path.join(root, f)
            print("FILE:", os.path.relpath(path, OUT_DIR), os.path.getsize(path))


if __name__ == "__main__":
    main()

