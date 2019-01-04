import numpy as np

def compute_fdr(df, args):

    fdr = args["false_discovery_rate"]

    df = df.sort_values("NLP", ascending=False)

    unique_nlp = df.NLP.drop_duplicates()

    sizes = unique_nlp.apply(lambda x: (df.NLP == x).sum())
    indices = unique_nlp.apply(lambda x: (df.NLP >= x).sum())

    nlrs = []
    for nlp, j in zip(unique_nlp, indices):
        m = sum(df.MAX_NLP >= nlp)
        by = np.log10(1/sum(range(1, m + 1)))
        nls = nlp + np.log10(j/m)
        nlrs.append(max(nls - by, 0))


    nlqs = [max(nlrs[i:]) for i in range(len(nlrs))]

    nlqss = np.repeat(nlqs, sizes)

    qvals = 10 ** -nlqss

    df.insert(df.shape[1], "QVAL", qvals)

    df = df[df.QVAL <= fdr]

    return df
