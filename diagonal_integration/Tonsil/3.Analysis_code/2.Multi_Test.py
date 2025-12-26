import pandas as pd
import numpy as np
import statsmodels.stats.multitest
import pyreadr

res = pyreadr.read_r('./multi_test.rds')[None]
intersections = pd.DataFrame.from_records(
    res["Intersections"].str.split(" & ").map(lambda x: {
        item.split("_")[1]: item.split("_")[0]
        for item in x
    }), index=res.index
)
res = res.join(intersections).dropna(
    subset=intersections.columns
)  # Dropped are those involving intersection within the same domains
res["FDR"] = statsmodels.stats.multitest.fdrcorrection(res["P.value"])[1]
res["-log10 P.value"] = -np.log10(res["P.value"])
res["-log10 FDR"] = -np.log10(res["FDR"])
res = res.sort_values("FDR")
res_samect = res.query("`scRNA` == `scRNAp` == `codex` == `codexp`")
res_samect.to_csv("./marker_intersection_test.csv",
           index=False)