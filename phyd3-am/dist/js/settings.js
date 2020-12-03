const summary_table_json = "data/summary_table.json"
const validation_xml = "data/SARS-CoV-2.xml"
const metadata_json = "data/SARS-CoV-2.xml.json"
const stats_json = "data/SARS-CoV-2.xml.stats.json"
const assay_result_path = "data/assay_result_json"
const db_totals_json = "data/db_totals.json"
const country_latlngs_json = "country_latlngs.json"
const assay_val_results_json = "data/SARS-CoV-2.xml.geo.json"

const phyd3_opts = {
    dynamicHide: true,
    // height: 2400,
    invertColors: false,
    lineupNodes: true,
    showDomains: true,
    showDomainNames: false,
    showDomainColors: true,
    showGraphs: true,
    showGraphLegend: true,
    showLength: false,
    showNodeNames: true,
    showNodesType: "only leaf",
    showPhylogram: false,
    showTaxonomy: true,
    showFullTaxonomy: false,
    showSequences: false,
    showTaxonomyColors: true,
    backgroundColor: "#FFFFFF",
    foregroundColor: "#000000",
    nanColor: "#f5f5f5"
};

const metadata = $.getJSON(metadata_json)
const db_stats = $.getJSON(db_totals_json)
const assay_stats = $.getJSON(stats_json)
const summary_json = $.getJSON(summary_table_json)
const country_latlngs = $.getJSON(country_latlngs_json)
const assay_val_results = $.getJSON(assay_val_results_json)
