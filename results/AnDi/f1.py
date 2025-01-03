#!/usr/bin/env python3
import numpy as np
import pandas as pd
from sklearn.metrics import (f1_score, multilabel_confusion_matrix)

# multilabel_confusion_matrix(y_true, y_pred)

# identical to mean accuracy
def F1(y_true, y_pred):
    return f1_score(y_true, y_pred, average="micro")

for fname in ["pred", "pred_134"]:
    print(fname)
    df = pd.read_table(f"./matroid_generators/CNN/{fname}.tsv")
    f1 = F1(df.Label, df.Pred)
    np.mean(df.Label == df.Pred)
    print("F1 =", f1)
    bootstraps = []
    for i in range(100):
        y_pred = np.random.choice(np.unique(df.Label), len(df.Label))
        f1 = F1(df.Label, y_pred)
        bootstraps.append(f1)
    print("random guess F1 =", np.mean(bootstraps))

