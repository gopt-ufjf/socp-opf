import pandas as pd
import pyomo.environ as pyo
import numpy as np
from collections import defaultdict

import pandas as pd
import numpy as np
import pyomo.environ as pyo

def create_report(model):

    vars_1d = {}
    vars_2d = {}
    vars_3d = {}

    # percorre todas as variáveis do modelo
    for var in model.component_objects(pyo.Var, active=True):
        v = getattr(model, var.name)

        for idx in v:
            val = pyo.value(v[idx])

            # transforma índice em tupla
            if not isinstance(idx, tuple):
                idx = (idx,)

            dim = len(idx)

            if dim == 1:
                vars_1d.setdefault(var.name, {})[idx[0]] = val

            elif dim == 2:
                vars_2d.setdefault(var.name, {})[idx] = val

            elif dim == 3:
                vars_3d.setdefault(var.name, {})[idx] = val

    # ---- DataFrame 1 índice ----
    df1 = None
    if vars_1d:
        all_index = sorted(
            set(i for d in vars_1d.values() for i in d.keys())
        )

        df1 = pd.DataFrame(index=all_index)

        for var, data in vars_1d.items():
            df1[var] = pd.Series(data)

    # ---- DataFrame 2 índices ----
    df2 = None
    if vars_2d:
        all_index = sorted(
            set(i for d in vars_2d.values() for i in d.keys())
        )

        df2 = pd.DataFrame(index=pd.MultiIndex.from_tuples(all_index))

        for var, data in vars_2d.items():
            df2[var] = pd.Series(data)

    # ---- DataFrame 3 índices ----
    df3 = None
    if vars_3d:
        all_index = sorted(
            set(i for d in vars_3d.values() for i in d.keys())
        )

        df3 = pd.DataFrame(index=pd.MultiIndex.from_tuples(all_index))

        for var, data in vars_3d.items():
            df3[var] = pd.Series(data)
    else:
        df3 = pd.DataFrame()

    return df1, df2, df3



