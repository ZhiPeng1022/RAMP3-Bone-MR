# 解压破骨细胞 GEO 数据集并列出文件

import os
import tarfile


DATA = r"data/derived"
TARS = ["GSE303003_RAW.tar", "GSE310113_RAW.tar"]


def main():
    for tar_name in TARS:
        tar_path = os.path.join(DATA, tar_name)
        out_dir = os.path.join(DATA, tar_name.replace(".tar", "_unpacked"))
        if not os.path.exists(tar_path):
            print("缺少文件：", tar_path)
            continue
        os.makedirs(out_dir, exist_ok=True)
        with tarfile.open(tar_path) as tf:
            members = tf.getmembers()
            print("=====", tar_name, "members:", len(members))
            for m in members[:40]:
                print(" ", m.name, m.size)
            tf.extractall(out_dir)
        print("解压到：", out_dir)


if __name__ == "__main__":
    main()

