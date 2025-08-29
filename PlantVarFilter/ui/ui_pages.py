# ui/ui_pages.py
import dearpygui.dearpygui as dpg


def page_preprocess_samtools(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_pre_sam"):
        dpg.add_text("\nClean BAM: sort / fixmate / markdup / index + QC reports", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                bam_btn = dpg.add_button(
                    label="Choose a BAM file",
                    callback=lambda: dpg.show_item("file_dialog_bam"),
                    width=220,
                    tag="tooltip_bam_sam",
                )
                app._secondary_buttons.append(bam_btn)
                dpg.add_text("No file", tag="sam_bam_path_lbl", wrap=500)

            with dpg.group():
                app.sam_threads = dpg.add_input_int(
                    label="Threads", width=220, default_value=4, min_value=1, min_clamped=True
                )
                app._inputs.append(app.sam_threads)

                app.sam_remove_dups = dpg.add_checkbox(
                    label="Remove duplicates (instead of marking)", default_value=False
                )
                app._inputs.append(app.sam_remove_dups)

                app.sam_compute_stats = dpg.add_checkbox(
                    label="Compute QC reports (flagstat/stats/idxstats/depth)", default_value=True
                )
                app._inputs.append(app.sam_compute_stats)

                dpg.add_spacer(height=6)
                app.sam_out_prefix = dpg.add_input_text(
                    label="Output prefix (optional)",
                    hint="Leave empty to auto-generate next to the input file",
                    width=320,
                )
                app._inputs.append(app.sam_out_prefix)

                dpg.add_spacer(height=12)
                run_sam = dpg.add_button(
                    label="Run samtools preprocess",
                    callback=app.run_samtools_preprocess,
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_sam)
    return "page_pre_sam"


def page_variant_calling(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_vc"):
        dpg.add_text("\nCall variants with bcftools mpileup + call", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                b1 = dpg.add_button(
                    label="Choose BAM (single)",
                    callback=lambda: dpg.show_item("file_dialog_bam_vc"),
                    width=220,
                    tag="tooltip_bam_vc",
                )
                app._secondary_buttons.append(b1)
                dpg.add_text("", tag="vc_bam_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                b2 = dpg.add_button(
                    label="Choose BAM-list (.list)",
                    callback=lambda: dpg.show_item("file_dialog_bamlist"),
                    width=220,
                    tag="tooltip_bamlist_vc",
                )
                app._secondary_buttons.append(b2)
                dpg.add_text("", tag="vc_bamlist_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                fa_btn = dpg.add_button(
                    label="Choose reference FASTA",
                    callback=lambda: dpg.show_item("file_dialog_fasta"),
                    width=220,
                    tag="tooltip_fa_vc",
                )
                app._secondary_buttons.append(fa_btn)
                dpg.add_text("", tag="vc_ref_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                reg_btn2 = dpg.add_button(
                    label="Choose regions BED (optional)",
                    callback=lambda: dpg.show_item("file_dialog_blacklist"),
                    width=220,
                    tag="tooltip_reg_vc",
                )
                app._secondary_buttons.append(reg_btn2)
                dpg.add_text("", tag="vc_regions_path_lbl", wrap=500)

            with dpg.group():
                app.vc_threads = dpg.add_input_int(
                    label="Threads", width=220, default_value=4, min_value=1, min_clamped=True
                )
                app._inputs.append(app.vc_threads)

                app.vc_ploidy = dpg.add_input_int(
                    label="Ploidy",
                    width=220,
                    default_value=2,
                    min_value=1,
                    min_clamped=True,
                    tag="tooltip_ploidy",
                )
                app._inputs.append(app.vc_ploidy)

                app.vc_min_bq = dpg.add_input_int(
                    label="Min BaseQ",
                    width=220,
                    default_value=20,
                    min_value=0,
                    min_clamped=True,
                    tag="tooltip_bq",
                )
                app._inputs.append(app.vc_min_bq)

                app.vc_min_mq = dpg.add_input_int(
                    label="Min MapQ",
                    width=220,
                    default_value=20,
                    min_value=0,
                    min_clamped=True,
                    tag="tooltip_mq",
                )
                app._inputs.append(app.vc_min_mq)

                dpg.add_spacer(height=6)
                app.vc_out_prefix = dpg.add_input_text(
                    label="Output prefix (optional)",
                    hint="Leave empty to auto-generate next to the BAM",
                    width=320,
                )
                app._inputs.append(app.vc_out_prefix)

                dpg.add_spacer(height=12)
                run_vc = dpg.add_button(
                    label="Call variants (bcftools)",
                    callback=app.run_variant_calling,
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_vc)
    return "page_vc"


def page_preprocess_bcftools(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_pre_bcf"):
        dpg.add_text("\nNormalize / split multiallelic / sort / filter / set IDs (bcftools)", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                vcf_btn_bcf = dpg.add_button(
                    label="Choose a VCF file",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                    tag="tooltip_vcf_bcf",
                )
                app._secondary_buttons.append(vcf_btn_bcf)
                dpg.add_text("", tag="bcf_vcf_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                fasta_btn = dpg.add_button(
                    label="Choose reference FASTA (for left-align)",
                    callback=lambda: dpg.show_item("file_dialog_fasta"),
                    width=220,
                    tag="tooltip_fa_bcf",
                )
                app._secondary_buttons.append(fasta_btn)
                dpg.add_text("", tag="bcf_ref_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                reg_btn = dpg.add_button(
                    label="Choose regions BED (optional)",
                    callback=lambda: dpg.show_item("file_dialog_blacklist"),
                    width=220,
                    tag="tooltip_reg_bcf",
                )
                app._secondary_buttons.append(reg_btn)
                dpg.add_text("", tag="bcf_regions_path_lbl", wrap=500)

            with dpg.group():
                app.bcf_split = dpg.add_checkbox(label="Split multiallelic", default_value=True)
                app._inputs.append(app.bcf_split)

                app.bcf_left = dpg.add_checkbox(label="Left-align indels (needs FASTA)", default_value=True)
                app._inputs.append(app.bcf_left)

                app.bcf_sort = dpg.add_checkbox(label="Sort", default_value=True)
                app._inputs.append(app.bcf_sort)

                app.bcf_setid = dpg.add_checkbox(label="Set ID to CHR:POS:REF:ALT", default_value=True)
                app._inputs.append(app.bcf_setid)

                app.bcf_compr = dpg.add_checkbox(label="Compress output (.vcf.gz)", default_value=True)
                app._inputs.append(app.bcf_compr)

                app.bcf_index = dpg.add_checkbox(label="Index output (tabix)", default_value=True)
                app._inputs.append(app.bcf_index)

                app.bcf_rmflt = dpg.add_checkbox(label="Keep only PASS (remove filtered)", default_value=False)
                app._inputs.append(app.bcf_rmflt)

                dpg.add_spacer(height=6)
                app.bcf_filter_expr = dpg.add_input_text(
                    label="bcftools filter expression (optional)",
                    hint="Example: QUAL>=30 && INFO/DP>=10",
                    width=320,
                )
                app._inputs.append(app.bcf_filter_expr)

                dpg.add_spacer(height=6)
                app.bcf_out_prefix = dpg.add_input_text(
                    label="Output prefix (optional)",
                    hint="Leave empty to auto-generate next to the input file",
                    width=320,
                )
                app._inputs.append(app.bcf_out_prefix)

                dpg.add_spacer(height=12)
                run_bcf = dpg.add_button(
                    label="Run bcftools preprocess",
                    callback=app.run_bcftools_preprocess,
                    width=240,
                    height=38,
                )
                app._primary_buttons.append(run_bcf)
    return "page_pre_bcf"


def page_check_vcf(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_check_vcf"):
        dpg.add_text("\nCheck VCF quality before conversion/analysis", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                vcf_btn_qc = dpg.add_button(
                    label="Choose a VCF file",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                    tag="tooltip_vcf_qc",
                )
                app._secondary_buttons.append(vcf_btn_qc)
                dpg.add_text("", tag="qc_vcf_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                bl_btn = dpg.add_button(
                    label="Choose blacklist BED (optional)",
                    callback=lambda: dpg.show_item("file_dialog_blacklist"),
                    width=220,
                    tag="tooltip_bl_qc",
                )
                app._secondary_buttons.append(bl_btn)
                dpg.add_text("", tag="qc_bl_path_lbl", wrap=500)

            with dpg.group():
                app.deep_scan = dpg.add_checkbox(label="Deep scan", default_value=False)
                app._inputs.append(app.deep_scan)

                dpg.add_spacer(height=12)
                run_qc_btn = dpg.add_button(
                    label="Run Quality Check",
                    callback=app.run_vcf_qc,
                    width=200,
                    height=36,
                )
                app._primary_buttons.append(run_qc_btn)
    return "page_check_vcf"


def page_convert_plink(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_plink"):
        dpg.add_text("\nConvert a VCF file into PLINK BED and apply MAF/missing genotype filters.", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                dpg.add_text("Select files:", indent=0)
                vcf = dpg.add_button(
                    label="Choose a VCF file",
                    callback=lambda: dpg.show_item("file_dialog_vcf"),
                    width=220,
                    tag="tooltip_vcf",
                )
                app._secondary_buttons.append(vcf)
                dpg.add_text("", tag="conv_vcf_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                variant_ids = dpg.add_button(
                    label="Choose IDs file (option)",
                    callback=lambda: dpg.show_item("file_dialog_variants"),
                    width=220,
                    tag="tooltip_variant",
                )
                app._secondary_buttons.append(variant_ids)
                dpg.add_text("", tag="ids_path_lbl", wrap=500)

            with dpg.group():
                dpg.add_text("Apply filters:", indent=0)
                maf_input = dpg.add_input_float(
                    label="Minor allele frequency (MAF)",
                    width=220,
                    default_value=0.05,
                    step=0.005,
                    tag="tooltip_maf",
                )
                app._inputs.append(maf_input)

                dpg.add_spacer(height=6)
                geno_input = dpg.add_input_float(
                    label="Missing genotype rate",
                    width=220,
                    default_value=0.10,
                    step=0.005,
                    tag="tooltip_missing",
                )
                app._inputs.append(geno_input)

                dpg.add_spacer(height=14)
                convert_btn = dpg.add_button(
                    label="Convert VCF",
                    callback=lambda s, a, u=[maf_input, geno_input, vcf, variant_ids]: app.convert_vcf(s, a, u),
                    width=160,
                    height=36,
                )
                app._primary_buttons.append(convert_btn)
    return "page_plink"


def page_gwas(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_gwas"):
        dpg.add_text("\nStart GWAS Analysis", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                geno = dpg.add_button(
                    label="Choose a BED file",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=220,
                    tag="tooltip_bed",
                )
                app._secondary_buttons.append(geno)
                dpg.add_text("", tag="gwas_bed_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                pheno = dpg.add_button(
                    label="Choose a phenotype file",
                    callback=lambda: dpg.show_item("file_dialog_pheno"),
                    width=220,
                    tag="tooltip_pheno",
                )
                app._secondary_buttons.append(pheno)
                dpg.add_text("", tag="gwas_pheno_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                cov_file = dpg.add_button(
                    label="Choose covariate file (option)",
                    callback=lambda: dpg.show_item("file_dialog_cov"),
                    width=220,
                    tag="tooltip_cov",
                )
                app._secondary_buttons.append(cov_file)
                dpg.add_text("", tag="gwas_cov_path_lbl", wrap=500)

            with dpg.group():
                app.gwas_combo = dpg.add_combo(
                    label="Analysis Algorithms",
                    items=[
                        "FaST-LMM",
                        "Linear regression",
                        "Ridge Regression",
                        "Random Forest (AI)",
                        "XGBoost (AI)",
                    ],
                    width=240,
                    default_value="FaST-LMM",
                    tag="tooltip_algorithm",
                )
                app._inputs.append(app.gwas_combo)

                dpg.add_spacer(height=14)
                gwas_btn = dpg.add_button(
                    label="Run GWAS",
                    callback=lambda s, a, u=[geno, pheno]: app.run_gwas(s, a, u),
                    width=180,
                    height=36,
                )
                app._primary_buttons.append(gwas_btn)
    return "page_gwas"


def page_genomic_prediction(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_gp"):
        dpg.add_text("\nStart Genomic Prediction", indent=10)
        dpg.add_spacer(height=10)
        with dpg.group(horizontal=True, horizontal_spacing=60):
            with dpg.group():
                geno = dpg.add_button(
                    label="Choose a BED file",
                    callback=lambda: dpg.show_item("file_dialog_bed"),
                    width=220,
                    tag="tooltip_bed_gp",
                )
                app._secondary_buttons.append(geno)
                dpg.add_text("", tag="gp_bed_path_lbl", wrap=500)

                dpg.add_spacer(height=6)
                pheno = dpg.add_button(
                    label="Choose a phenotype file",
                    callback=lambda: dpg.show_item("file_dialog_pheno"),
                    width=220,
                    tag="tooltip_pheno_gp",
                )
                app._secondary_buttons.append(pheno)
                dpg.add_text("", tag="gp_pheno_path_lbl", wrap=500)

            with dpg.group():
                app.gwas_gp = dpg.add_combo(
                    label="Analysis Algorithms",
                    items=[
                        "XGBoost (AI)",
                        "Random Forest (AI)",
                        "Ridge Regression",
                        "GP_LMM",
                        "val",
                    ],
                    width=240,
                    default_value="XGBoost (AI)",
                    tag="tooltip_algorithm_gp",
                )
                app._inputs.append(app.gwas_gp)

                dpg.add_spacer(height=14)
                gp_btn = dpg.add_button(
                    label="Run Genomic Prediction",
                    callback=lambda s, a, u=[geno, pheno]: app.run_genomic_prediction(s, a, u),
                    width=220,
                    height=36,
                )
                app._primary_buttons.append(gp_btn)
    return "page_gp"


def page_settings(app, parent):
    with dpg.group(parent=parent, show=False, tag="page_settings"):
        dpg.add_spacer(height=10)

        with dpg.table(
            header_row=False,
            borders_innerH=False,
            borders_outerH=False,
            borders_innerV=False,
            borders_outerV=False,
            resizable=False,
        ):
            dpg.add_table_column()
            dpg.add_table_column()

            with dpg.table_row():
                with dpg.group():
                    dpg.add_text("General Settings", color=(200, 180, 90))
                    dpg.add_spacer(height=8)

                    dpg.add_checkbox(
                        label="Night Mode (Dark)",
                        default_value=True,
                        tag="settings_dark_toggle",
                    )

                    dpg.add_spacer(height=6)
                    app.nr_jobs = dpg.add_input_int(
                        label="Number of jobs to run",
                        width=220,
                        default_value=-1,
                        step=1,
                        min_value=-1,
                        max_value=50,
                        min_clamped=True,
                        max_clamped=True,
                        tag="tooltip_nr_jobs",
                    )
                    app._inputs.append(app.nr_jobs)

                    app.gb_goal = dpg.add_input_int(
                        label="Gigabytes of memory per run",
                        width=220,
                        default_value=0,
                        step=4,
                        min_value=0,
                        max_value=512,
                        min_clamped=True,
                        max_clamped=True,
                        tag="tooltip_gb_goal",
                    )
                    app._inputs.append(app.gb_goal)

                    app.plot_stats = dpg.add_checkbox(
                        label="Advanced Plotting", default_value=False, tag="tooltip_stats"
                    )
                    app._inputs.append(app.plot_stats)

                    app.snp_limit = dpg.add_input_text(
                        label="SNP limit", width=220, default_value="", tag="tooltip_limit"
                    )
                    app._inputs.append(app.snp_limit)

                with dpg.group():
                    dpg.add_text("Machine Learning Settings", color=(200, 180, 90))
                    dpg.add_spacer(height=8)

                    app.train_size_set = dpg.add_input_int(
                        label="Training size %",
                        width=220,
                        default_value=70,
                        step=10,
                        min_value=0,
                        max_value=100,
                        min_clamped=True,
                        max_clamped=True,
                        tag="tooltip_training",
                    )
                    app._inputs.append(app.train_size_set)

                    app.estim_set = dpg.add_input_int(
                        label="Number of trees",
                        width=220,
                        default_value=200,
                        step=10,
                        min_value=1,
                        min_clamped=True,
                        tag="tooltip_trees",
                    )
                    app._inputs.append(app.estim_set)

                    app.max_dep_set = dpg.add_input_int(
                        label="Max depth",
                        width=220,
                        default_value=3,
                        step=10,
                        min_value=0,
                        max_value=100,
                        min_clamped=True,
                        max_clamped=True,
                        tag="tooltip_depth",
                    )
                    app._inputs.append(app.max_dep_set)

                    app.model_nr = dpg.add_input_int(
                        label="Nr. of models",
                        width=220,
                        default_value=1,
                        step=1,
                        min_value=1,
                        max_value=50,
                        min_clamped=True,
                        max_clamped=True,
                        tag="tooltip_model",
                    )
                    app._inputs.append(app.model_nr)

                    app.aggregation_method = dpg.add_combo(
                        ("sum", "median", "mean"),
                        label="Aggregation Method",
                        width=220,
                        default_value="sum",
                        tag="tooltip_aggr",
                    )
                    app._inputs.append(app.aggregation_method)

        dpg.add_spacer(height=18)
        dpg.add_text("Appearance", color=(200, 180, 90))
        dpg.add_spacer(height=8)

        dpg.add_slider_float(
            label="Font scale",
            min_value=0.85,
            max_value=1.35,
            default_value=1.10,
            width=520,
            tag="settings_font_scale",
        )

        dpg.add_spacer(height=8)

        dpg.add_combo(
            items=[
                "Evergreen (Green)",
                "Teal",
                "Blue",
                "Amber",
                "Purple",
            ],
            default_value="Evergreen (Green)",
            width=260,
            label="Accent color",
            tag="settings_accent_combo",
        )
    return "page_settings"


def build_pages(app, parent):
    pages = {}
    pages["pre_sam"] = page_preprocess_samtools(app, parent)
    pages["vc"] = page_variant_calling(app, parent)
    pages["pre_bcf"] = page_preprocess_bcftools(app, parent)
    pages["check_vcf"] = page_check_vcf(app, parent)
    pages["plink"] = page_convert_plink(app, parent)
    pages["gwas"] = page_gwas(app, parent)
    pages["gp"] = page_genomic_prediction(app, parent)
    pages["settings"] = page_settings(app, parent)
    return pages
