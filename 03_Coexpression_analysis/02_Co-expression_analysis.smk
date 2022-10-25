# Gene co-expression analysis by Bingru Zhao    #

# 1. Method: WGCNA
rule run_wgcna:
    script:
        "WGCNA.R"

# 2. Method: MEGENA
rule run_megena:
    script:
        "MEGENA.R"

# 3. Method: ICA
rule run_ica:
    script:
        "ICA.R"

# 4. Method: CEMiTool
rule run_cemitool:
    script:
        "CEMiTool.R"

# 5. Method: PEER
rule run_peer:
    script:
        "PEER.R"

