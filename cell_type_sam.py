import pandas as pd

folder = "/path/to/folder/with/cell/barcode/csvs/"
files = ("B-cells", "Malignant", "DC", "Macro", "Mono", "NK", "Other", "pDC",
         "Plasma", "Th", "Tcyt", "Treg") # the changing tags in the names of the csv files

for file in files:
    print(file)
    df = pd.read_csv(folder+file+".csv", header = 1, names = ["drop", "barcode"])
    df = df.drop(["drop"], axis=1)
    df.barcode = df.barcode.str.replace('P12_1670_', '')
    df.barcode = df.barcode.str.replace('-1', '')
    barcodelist = df.barcode.tolist()
    outfile = "/output/path/{}.sam".format(file)
    infile = "/path/to/samfile/samfile.sam"
    with open(outfile, "a") as outsam:
        with open(infile) as insam:
            for line in insam:
                if not line.startswith("@"):
                    a = line
                    line = line.split("\t")
                    d = {}
                    for elem in line:
                        if len(elem) > 7:
                            d[elem[:5]] = elem[5:-2]
                    if "CB:Z:" in d:
                        if d["CB:Z:"] in barcodelist:
                            outsam.write(a)
                else:
                    outsam.write(line)
print("Finished")