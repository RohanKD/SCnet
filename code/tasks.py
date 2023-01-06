from tabulate import tabulate

import anndata

if __name__ == '__main__':
    adata = anndata.read_h5ad("../data/Mouse_pFC.h5ad")
    # samples = []
    # table = []
    # table2 = []
    # print(adata)
    # sample = str(adata.obs['Experiment']).split()
    # print(adata.n_vars)
    data = adata[adata.obs['Sample'] == "PFCSample11"]
    print(data.obs)
    print(data.var)
    # table.append(str(data.obs['cell.type'].value_counts()).split()[0:16])
    # data = adata[adata.obs['Experiment'] == "Cortex2"]
    # table.append(str(data.obs['cell.type'].value_counts()).split()[0:14])
    # print(tabulate(table, tablefmt="fancy_grid"))

# for i in sample:
#     if i.count("P") >= 1:
#         if i.count("_") == 0:
#             if 7 < len(i) < 20:
#                 if samples.count(i.replace("[", "").replace("]", "").replace(",", "").replace("'", "")) == 0:
#                     samples.append(i.replace("[", "").replace("]", "").replace(",", "").replace("'", ""))
# header =["Sample", "Cell", "#", "Cell", "#", "Cell", "#","Cell", "#","Cell", "#","Cell", "#","Cell", "#"]
# #table.append(header)
# for j in samples:
#     data = adata[adata.obs['sampleID'] == j]
#     row = str(data.obs['cell.type'].value_counts()).split()[0:18]
#     print(j)
#     row.insert(0, j)
#     print(row)
#     table.append(row)
#
# print(tabulate(table, tablefmt="fancy_grid"))

# counter = 0
# for i in sample:
#     if i.count("sampleID") == 1:
#         if i.count("_") == 0:
#             if 7 < len(i) < 20:
#                 if samples.count(i.replace("[", "").replace("]", "").replace(",", "").replace("'", "")) == 0:
#                     samples.append(i.replace("[", "").replace("]", "").replace(",", "").replace("'", ""))
# header = ["Sample", "Cell", "Count", "Cell", "Count", "Cell", "Count", "Cell", "Count", "Cell", "Count", "Cell",
#           "Count", "Cell", "Count", "Total"]
# table.append(header)
# # print(samples)
# for j in samples:
#     counter += 1
#
#     # if (counter%2 == 0):
#     # print(j)
#     sample1_adata = adata[adata.obs['sampleID'] == j]
#     row = str(sample1_adata.obs['cell.type'].value_counts()).split()
#     row = row[0:14]
#     row_b = row
#     total = 0
#     counter = 0
#     for r in row:
#         counter += 1
#         if counter % 2 == 0:
#             # print(r)
#             total += int(r)
#     counter = 1
#     for p in range(7):
#         #row_b[counter] = str(round(int(row[counter]) / total, 4))
#         counter += 2
#
#     row.insert(0, j)
#     row.insert(len(row), str(total))
#     table.append(row)
#     table2.append(row_b)
#
# print(tabulate(table, tablefmt="fancy_grid"))
# # print(tabulate(table2, tablefmt="fancy_grid"))
#
# # print(sample)
#
# # print(adata[adata.obs['Sample'] == "PFCSample9"].obs['cell.type'].value_counts())
#
# total_cells = []
# props = []
