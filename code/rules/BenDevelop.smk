rule CollectSummaryStatsForColoc:
    input:
        expand("hyprcoloc/LociWiseSummaryStatsInput/{FeatureCoordinatesRedefinedFor}", FeatureCoordinatesRedefinedFor=["ForColoc", "ForGWASColoc"])

