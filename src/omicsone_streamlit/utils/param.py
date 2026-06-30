import numpy as np
import pandas as pd

import os, sys, re
import configparser


class Parameter:
    def __init__(self, path) -> None:
        self.path = path
        self.read_param(path)

    def __getitem__(self, key):
        if key not in self.__dict__.keys():
            raise KeyError
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def read_param(self, path):
        with open(path, "r") as f:
            lines = f.readlines()
            for l in lines:
                l = str(l).strip()
                if re.search("^#", l):
                    continue
                if l == "":
                    continue
                key, value = re.split("=", l)
                self[key] = value


class DiffParams:
    def __init__(self):
        # Initialize DiffSettings attributes to None by default
        self.sources = None
        self.target_column = None
        self.labels = None
        self.missing_process = None
        self.missing_samples_ratio = None
        self.method = None
        self.fdr = None
        self.fold_change = None

    def read_from_config(self, config):
        # # get source matrix
        self.sources = config.get("DiffSettings", "Sources", fallback="All")
        self.target_column = config.get(
            "DiffSettings", "TargetColumn", fallback="SampleStatus"
        )
        self.labels = config.get("DiffSettings", "Labels", fallback="T,N").split(",")
        self.missing_process = config.get(
            "DiffSettings", "MissingProcess", fallback="Remove"
        )
        self.missing_samples_ratio = config.getint(
            "DiffSettings", "MissingSamplesRatio", fallback="50"
        )
        self.method = config.get(
            "DiffSettings", "Method", fallback="Wilcoxon(Unpaired)"
        )
        self.fdr = config.getfloat("DiffSettings", "FDR", fallback="0.01")
        self.fold_change = config.getint("DiffSettings", "FoldChange", fallback="2")


class ProfileParams:
    def __init__(self):
        self.target_column = None
        self.corr_method = None

    def read_from_config(self, config):
        self.target_column = config.get(
            "ProfileSettigs", "TargetColumn", fallback="SampleStatus"
        )
        self.corr_method = config.get(
            "ProfileSettings", "CorrMethod", fallback="Spearman"
        )


class PreprocessParams:
    def __init__(self):
        self.normalization_method = None

    def read_from_config(self, config):
        self.normalization_method = config.get(
            "PreprocessSettings", "NormalizationMethod", fallback="MedianNorm"
        )


class DecompParams:
    def __init__(self):
        self.target_column = None

    def read_from_config(self, config):
        self.target_column = config.get(
            "DecompositionSettings", "TargetColumn", fallback="SampleStatus"
        )


class ClusterParams:
    def __init__(self):
        self.target_column = None

    def read_from_config(self, config):
        self.target_column = config.get(
            "ClusteringSettings", "TargetColumn", fallback="SampleStatus"
        )


class EnrichParams:
    def __init__(self):
        self.reference = None

    def read_from_config(self, config):
        self.reference = config.get(
            "EnrichmentSettings", "Reference", fallback="KEGG2016"
        )


class SurvivalParams:
    def __init__(self):
        self.sources = None
        self.vital_column = None
        self.days_column = None

    def read_from_config(self, config):
        self.sources = config.get("SurvivalSettings", "Sources", fallback="All")
        self.vital_column = config.get(
            "SurvivalSettings", "VitalColumn", fallback="VitalStatus"
        )
        self.days_column = config.get(
            "SurvivalSettings", "SurvivalDaysColumn", fallback="SurvivalDays"
        )


class GlycosylationParams:
    def __init__(self):
        self.genes = None

    def read_from_config(self, config):
        self.genes = config.get("GlycosylationSettings", "Genes", fallback="").split(
            "\n"
        )
        self.genes = [gene.strip() for gene in self.genes]


class CrossTalkParams:
    def __init__(self):
        self.ptm1 = None
        self.ptm2 = None
        self.gsea_ora_gene_sets = None

    def read_from_config(self, config):
        self.ptm1 = config.get("CrossTalkSettings", "PTM1", fallback="Glycosylation")
        self.ptm2 = config.get("CrossTalkSettings", "PTM2", fallback="Phosphorylation")
        gene_sets = config.get("CrossTalkSettings", "GSEA_ORA_GeneSets", fallback="MSigDB_Hallmark_2020")
        self.gsea_ora_gene_sets = gene_sets.split(",")

    def get_gene_sets(self):
        return self.gsea_ora_gene_sets


class AssociationParams:
    def __init__(self):
        self.target_column = None

    def read_from_config(self, config):
        self.target_column = config.get(
            "AssociationSettings", "TargetColumn", fallback="SampleStatus"
        )


class CommandParams:
    def __init__(self):
        self.meta_path = None
        self.protein_path = None
        self.glycopeptide_path = None
        self.phosphopeptide_path = None
        self.paths_map = None

        self.profile_card = False
        self.preprocess_card = False
        self.qc_card = False
        self.annotation_card = False
        self.diff_card = False
        self.decomposition_card = False
        self.clustering_card = False
        self.association_card = False
        self.enrich_card = False
        self.survival_card = False
        self.glycosylation_card = False
        self.crosstalk_card = False

    def get_matrix_path(self, name):
        if str(name).lower() == "protein":
            return self.protein_path
        elif str(name).lower() == "glycopeptide":
            return self.glycopeptide_path
        elif str(name).lower() == "phosphopeptide":
            return self.phosphopeptide_path
        
    def get_paths_map(self):
        if self.paths_map is None:
            # generate the paths map
            self.paths_map = dict()
            self.paths_map["meta"] = self.meta_path
            self.paths_map["protein"] = self.protein_path
            self.paths_map["glycopeptide"] = self.glycopeptide_path
            self.paths_map["phosphopeptide"] = self.phosphopeptide_path
        return self.paths_map

    def read_from_config_path(self, config_path):
        config = configparser.ConfigParser()
        config.read(config_path)

        self.meta_path = config.get("ExpressionMatrixPaths", "Meta")
        self.protein_path = config.get("ExpressionMatrixPaths", "Protein")
        self.glycopeptide_path = config.get("ExpressionMatrixPaths", "Glycopeptide")
        self.phosphopeptide_path = config.get("ExpressionMatrixPaths", "Phosphopeptide")

        self.profile_card = config.getboolean("Cards", "profile_card")
        self.preprocess_card = config.getboolean("Cards", "preprocess_card")
        self.qc_card = config.getboolean("Cards", "qc_card")
        self.annotation_card = config.getboolean("Cards", "annotation_card")
        self.diff_card = config.getboolean("Cards", "diff_card")
        self.decomposition_card = config.getboolean("Cards", "decomposition_card")
        self.clustering_card = config.getboolean("Cards", "clustering_card")
        self.association_card = config.getboolean("Cards", "association_card")
        self.enrich_card = config.getboolean("Cards", "enrich_card")
        self.survival_card = config.getboolean("Cards", "survival_card")
        self.glycosylation_card = config.getboolean("Cards", "glycosylation_card")
        self.crosstalk_card = config.getboolean("Cards", "crosstalk_card")

        self.diff_params = DiffParams()
        self.profile_params = ProfileParams()
        self.preprocess_params = PreprocessParams()
        self.decomp_params = DecompParams()
        self.cluster_params = ClusterParams()
        self.association_params = AssociationParams()
        self.enrich_params = EnrichParams()
        self.survival_params = SurvivalParams()
        self.gly_params = GlycosylationParams()
        self.crosstalk_params = CrossTalkParams()

        self.diff_params.read_from_config(config)
        self.profile_params.read_from_config(config)
        self.preprocess_params.read_from_config(config)
        self.decomp_params.read_from_config(config)
        self.cluster_params.read_from_config(config)
        self.association_params.read_from_config(config)
        self.enrich_params.read_from_config(config)
        self.survival_params.read_from_config(config)
        self.gly_params.read_from_config(config)
        self.crosstalk_params.read_from_config(config)
