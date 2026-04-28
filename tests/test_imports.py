import importlib

MODULES = [
    "sdi_pipeline",
    "sdi_robustness_checks",
    "effect_sign_enrichment",
    "gwas_pc_sensitivity",
    "candidate_gene_keyword_enrichment",
    "genotype_to_spectral_projection",
]


def test_modules_import():
    for module in MODULES:
        importlib.import_module(module)
