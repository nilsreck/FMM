
configurations = config.get("configurations", dict())

Control_list = configurations["Control_list"]
Normal_tissue_list = configurations["Normal_tissue_list"]
Cancer_list = configurations["Cancer_list"]
plot_bar = configurations["plot_bar"]
dimension_list = configurations["dimension_list"]
comparison_list = get_comparison_list(dimension_list)
max_iter = configurations["max_iter"]
annotation = configurations["Annotation"]
optimal_dim = configurations["optimal_dim"]

data_path = os.path.abspath("../Data")
resource_path = os.path.abspath("./resources")
#network_host = configurations["network_host"]
NCBI_account = configurations["NCBI_account"]
network_wrapper = configurations["network_access_wrapper"]

rule generate_networks:
    input:
        genes=get_genes(Control_list, Normal_tissue_list, Cancer_list, data_path),
        adj=data_path+"/Human_Biogrid_Adj_PPI_.npy",
        biogrid=data_path+"/Human_Biogrid_Genes_PPI_",
    output:
        networks=get_networks(Control_list, Normal_tissue_list, Cancer_list, data_path),
    params:
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        data_path=data_path,
    resources:
        mem_mb="8gb",
        nodes=1,
        runtime=1000,
    log:
        "logs/generate_networks.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/generate_networks.py"


rule generate_PPMI:
    input:
        networks=get_networks(Control_list, Normal_tissue_list, Cancer_list, data_path),
    output:
        matrices=get_PPMI_matrices(Control_list, Normal_tissue_list, Cancer_list, data_path),
    params:
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        data_path=data_path,
    resources:
        mem_mb="12gb",
        nodes=1,
        runtime=3000,
    log:
        "logs/generate_PPMI.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/generate_PPMI.py"


rule prepare_resources:
    output:
        dimensions=get_dimensions(dimension_list, resource_path),
        comparisons=get_comparisons(dimension_list, resource_path),
    params:
        dimension_list=dimension_list,
        resource_path=resource_path,
    resources:
        mem_mb="1gb",
        nodes=1,
        runtime=1000,
    log:
        "logs/prepare_resources.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/prepare_resources.py"


rule generate_embeddings:
    input:
        dimension=str(resource_path)+"/dimension_{dimension}",
        networks=get_networks(Control_list, Normal_tissue_list, Cancer_list, data_path),
        matrices=get_PPMI_matrices(Control_list, Normal_tissue_list, Cancer_list, data_path),
    output:
        canc_adj_embeddings=expand(data_path+"/Embeddings/_{matrix}_Matrix_{{dimension}}_PPI_{cancer}_Adj_Cancer.npy", cancer=Cancer_list, matrix=["P","U","G"]),
        canc_ppmi_embeddings=expand(data_path+"/Embeddings/_{matrix}_Matrix_{{dimension}}_PPI_{cancer}_PPMI_Cancer.npy", cancer=Cancer_list, matrix=["P","U","G"]),
        cont_adj_embeddings=expand(data_path+"/Embeddings/_{matrix}_Matrix_{{dimension}}_PPI_{control[1]}_{control[0]}_Adj_Control.npy", control=zip(Control_list, Normal_tissue_list), matrix=["P","U","G"]),
        cont_ppmi_embeddings=expand(data_path+"/Embeddings/_{matrix}_Matrix_{{dimension}}_PPI_{control[1]}_{control[0]}_PPMI_Control.npy", control=zip(Control_list, Normal_tissue_list), matrix=["P","U","G"]),
    params:
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        max_iter=max_iter,
        data_path=data_path,
    resources:
        mem_mb="18gb",
        nodes=4,
        runtime=20000,
    log:
        "logs/generate_embeddings{dimension}.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/generate_embeddings.py"


rule map_annotations:
    input:
        annotations=get_annotations(annotation, data_path),
        genes=get_genes(Control_list, Normal_tissue_list, Cancer_list, data_path),
        embeddings=get_embeddings(Control_list, Normal_tissue_list, Cancer_list, dimension_list, data_path),
    output:
        mapped_embeddings=get_mapped_embeddings(Control_list, Normal_tissue_list, Cancer_list, dimension_list, annotation, data_path),
    params:
        annotation=annotation,
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        dimension_list=dimension_list,
        data_path=data_path,
    resources:
        mem_mb="8gb",
        nodes=1,
        runtime=1000,
    log:
        "logs/map_annotations.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/map_annotations.py"


rule generate_FMMs:
    input:
        dimension=str(resource_path)+"/dimension_{dimension}",
        canc_adj_mapped_embeddings=expand(data_path+"/Embeddings/_GO_Embeddings_{annotation}_PPI_{cancer}_{{dimension}}_Adj_Cancer.csv", cancer=Cancer_list, annotation=annotation),
        canc_ppmi_mapped_embeddings=expand(data_path+"/Embeddings/_GO_Embeddings_{annotation}_PPI_{cancer}_{{dimension}}_PPMI_Cancer.csv", cancer=Cancer_list, annotation=annotation),
        cont_adj_mapped_embeddings=expand(data_path+"/Embeddings/_GO_Embeddings_{annotation}_PPI_{control[1]}_{control[0]}_{{dimension}}_Adj_Control.csv", control=zip(Control_list, Normal_tissue_list), annotation=annotation),
        cont_ppmi_mapped_embeddings=expand(data_path+"/Embeddings/_GO_Embeddings_{annotation}_PPI_{control[1]}_{control[0]}_{{dimension}}_PPMI_Control.csv", control=zip(Control_list, Normal_tissue_list), annotation=annotation),
    output:
        canc_fmms=expand(data_path+"/FMM/Cosine_Cancer_{cancer}_{{dimension}}_PPMI_{annotation}.csv", cancer=Cancer_list, annotation=annotation),
        cont_fmms=expand(data_path+"/FMM/Cosine_Control_{control[1]}_{control[0]}_{{dimension}}_PPMI_{annotation}.csv", control=zip(Control_list, Normal_tissue_list), annotation=annotation),
    params:
        annotation=annotation,
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        data_path=data_path,
    resources:
        mem_mb="12gb",
        nodes=min(6,len(Cancer_list)),
        runtime=1000,
    log:
        "logs/generate_FMMS{dimension}.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/generate_FMMs.py"


rule filter_FMMs:
    input:
        annotations=get_annotations(annotation, data_path),
        genes=get_genes(Control_list, Normal_tissue_list, Cancer_list, data_path),
    output:
        commons=get_common(Normal_tissue_list, annotation, data_path),
    params:
        annotation=annotation,
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        data_path=data_path,
    resources:
        mem_mb="6gb",
        nodes=1,
        runtime=2000,
    log:
        "logs/filter_FMMs.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/filter_FMMs.py"


rule functional_organization:
    input:
        fmms=get_FMMs(Control_list, Normal_tissue_list, Cancer_list, dimension_list, annotation, data_path),
        commons=get_common(Normal_tissue_list, annotation, data_path),
    output:
        get_similarity(Normal_tissue_list, annotation, data_path),
    params:
        annotation=annotation,
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        dimension_list=dimension_list,
        data_path=data_path,
    resources:
        mem_mb="2gb",
        nodes=1,
        runtime=1000,
    log:
        "logs/functional_organization.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/functional_organization.py"


rule optimal_dimensionality:
    input:
        #comparison=str(resource_path)+"/comparison_{comparison}",
        fmms=get_FMMs(Control_list, Normal_tissue_list, Cancer_list, dimension_list, annotation, data_path),
        bp=data_path+"/gene2go_Human_PPIGO_Specific_BP.json",
    output:
        canc_errors=expand("{data_path}/FMM/Relative_Cancer_{cancer}_PPMI_{annotation}.txt", data_path=data_path, cancer=Cancer_list, annotation=annotation),
        cont_errors=expand("{data_path}/FMM/Relative_Control_{control[1]}_{control[0]}_PPMI_{annotation}.txt", data_path=data_path, control=zip(Control_list, Normal_tissue_list), annotation=annotation),
    params:
        annotation=annotation,
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        dimension_list=dimension_list,
        data_path=data_path,
        #threads=6,
    resources:
        mem_mb="16gb",
        nodes=4,
        runtime=1000,
    log:
        "logs/optimal_dimensionality.log"
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/optimal_dimensionality.py"



rule plot_dimensionality_comparison:
    input:
        comparisons=get_comparisons(dimension_list, resource_path),
        canc_errors=expand("{data_path}/FMM/Relative_Cancer_{cancer}_PPMI_{annotation}.txt", data_path=data_path, cancer=Cancer_list, annotation=annotation, comparison=comparison_list),
        cont_errors=expand("{data_path}/FMM/Relative_Control_{control[1]}_{control[0]}_PPMI_{annotation}.txt", data_path=data_path, control=zip(Control_list, Normal_tissue_list), annotation=annotation, comparison=comparison_list),
    output:
        error_plot=get_error_plot(annotation, data_path),
    params:
        annotation=annotation,
        Cancer_list=Cancer_list,
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        dimension_list=dimension_list,
        data_path=data_path,
        #threads=6,
    resources:
        mem_mb="2gb",
        nodes=1,
        runtime=1000,
    log:
        "logs/plot_dimensionality_comparison.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/plot_dimensionality_comparison.py"

rule eval_movements:
    input:
        fmms=get_FMMs(Control_list, Normal_tissue_list, Cancer_list, dimension_list, annotation, data_path),
        commons=get_common(Normal_tissue_list, annotation, data_path),
        cancer_bp=data_path+"/enriched_in_cancer_gobp_terms_cosmic.txt",
    output:
        movements=get_movements(Normal_tissue_list, annotation, data_path),
        enrichments=data_path+"/cancer_enrichments_2std.svg",
    params:
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        optimal_dim=optimal_dim,
        data_path=data_path,
        annotation=annotation,
        plot_bar=plot_bar,
    resources:
        mem_mb="8gb",
        nodes=1,
        runtime=1000,
    log:
        "logs/eval_movements.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/eval_movements.py"


rule cancer_enrichment:
    input:
        cancer_bp=data_path+"/enriched_in_cancer_gobp_terms_cosmic.txt",
        file=expand("{data_path}/top_100_moving_{tissue}.csv", data_path=data_path,tissue=Normal_tissue_list),
        rank=expand("{data_path}/Rank_movement_{tissue}_PPMI_"+str(annotation)+".csv", data_path=data_path,tissue=Normal_tissue_list),
        common_genes=expand("{data_path}/Transformed_Common_{tissue}.csv", data_path=data_path, tissue=Normal_tissue_list)
    output:
        movement_tables=expand("{data_path}/Top_moving_{tissue}_Table.csv", data_path=data_path,tissue=Normal_tissue_list),
        counts=expand("{data_path}/Cancer_Count_{tissue}.csv", data_path=data_path, tissue=Normal_tissue_list),
    params:
        cmd=get_enrichment_cmd(network_wrapper, data_path, annotation, NCBI_account, Normal_tissue_list),
    resources:
        mem_mb="1gb",
        nodes=1,
        runtime=70000,
    log:
        get_enrichment_log()
    conda:
        "../../envs/enrichment.yaml"
    shell:
        """
        {params.cmd}
        """


rule analyze_functional_organization:
    input:
        genes=get_genes(Control_list, Normal_tissue_list, Cancer_list, data_path),
        embeddings=get_embeddings(Control_list, Normal_tissue_list, Cancer_list, dimension_list, data_path),
        mapped_embeddings=get_mapped_embeddings(Control_list, Normal_tissue_list, Cancer_list, dimension_list, annotation, data_path),
        movements=get_movements(Normal_tissue_list, annotation, data_path),
        annotations=get_annotations(annotation, data_path),
    output:
        distances=get_distances(Normal_tissue_list, data_path),
        organization=data_path+"/Organization_Common_Space.txt",
    params:
        Control_list=Control_list,
        Normal_tissue_list=Normal_tissue_list,
        Cancer_list=Cancer_list,
        optimal_dim=optimal_dim,
        data_path=data_path,
        annotation=annotation,
    resources:
        mem_mb="16gb",
        nodes=len(Cancer_list),
        runtime=20000,
    log:
        "logs/analyze_functional_organization.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/analyze_functional_organization.py"


rule gene_prediction:
    input:
        movements=get_movements(Normal_tissue_list, annotation, data_path),
        counts=get_counts(Normal_tissue_list, data_path),
        distances=get_distances(Normal_tissue_list, data_path),
    output:
        fold_rank=data_path+"/Fold_Rank_Table_Dos.txt",
        predictions=get_predictions(Normal_tissue_list, data_path),
    params:
        Normal_tissue_list=Normal_tissue_list,
        data_path=data_path,
    resources:
        mem_mb="2gb",
        nodes=1,
        runtime=1000,
    log:
        "logs/gene_prediction.log",
    conda:
        "../../envs/tools.yaml"
    script:
        "../scripts/gene_prediction.py"







