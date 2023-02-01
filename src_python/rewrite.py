import re
import sys

comment = "!"
single_index_kw = ("name", "from_file", "to_file", "transfer_func")
multi_index_kw = (
    "from_data_arrays",
    "target_coord_names",
    "upscale_ops",
    "limits",
)

# file names given as CLI args
file_in = sys.argv[1]
file_out = sys.argv[2]

in_block = False
block_no = 0

with open(file_in, "r") as f_in:
    lines = f_in.readlines()

with open(file_out, "w") as f_out:
    for line in lines:
        # copy comments
        if line.lstrip().startswith(comment):
            f_out.writelines([line])
            continue
        # copy anything outside data-array blocks
        if not in_block and not line.lstrip().startswith("name"):
            f_out.writelines([line])
            continue
        # check if we enter a block and if it is the first one
        if not in_block and line.lstrip().startswith("name"):
            in_block = True
            block_no += 1
        if in_block:
            # block ends with empty line
            if not line.strip():
                in_block = False
                f_out.writelines([line])
                continue
            # update single index
            if line.lstrip().startswith(single_index_kw):
                idx = line.index("=")
                sub = re.sub(r"\((\d+)\)", f"({block_no})", line[:idx])
                line = sub + line[idx:]
                f_out.writelines([line])
            # update latter multi-index
            elif line.lstrip().startswith(multi_index_kw):
                idx = line.index("=")
                sub = re.sub(r",(\d+)\)", f",{block_no})", line[:idx])
                line = sub + line[idx:]
                f_out.writelines([line])
            # write everything in between
            else:
                f_out.writelines([line])
