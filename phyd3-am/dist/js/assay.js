function updateStats(data, assay_stats){
    // Input JSON looks like this:
    // {
    //     "GISAID_tot": 13541,
    //     "GenBank_tot": 229,
    //     "final_date": "2020-04-27 00:00:00",
    //     "final_db_tot": 13558,
    //     "overlap": 212
    // }    
    $('#GISAID_tot').html(data.GISAID_tot.toLocaleString())
    $('#GenBank_tot').html(data.GenBank_tot.toLocaleString())
    $('#final_date').html(data.final_date)
    $('#overlap').html(data.overlap.toLocaleString())
    $('.final_db_tot').each(function(){
        $(this).html(data.final_db_tot.toLocaleString())
    })
    $('#final_col_tot').html(assay_stats.tree.collapsed_genome_num.toLocaleString())
    $('#final_leaf_tot').html(assay_stats.tree.leaf_num.toLocaleString())
}

function pyHeatmapAddPopper(){
    // Get user details for tooltip
    function nodeMetadata(){
        let id = $(this).html().split(" ").slice(-2)[0]
        let field_name = ""
        if( id.startsWith('EPI') ){
            field_name = "GISAID:"
        }
        html = `<div><div class="mt-2"><strong>${field_name} ${id}</strong></div><p>`

        $.each(metadata.responseJSON[id], function (idx, val) {
            if( val && val != "?" && val != "Unknown"){
                // strain  virus   gisaid_epi_isl  genbank_accession       date    
                // region  country division        location        
                // region_exposure country_exposure        division_exposure       segment length  
                // host    age     sex     pangolin_lineage        GISAID_clade    
                // originating_lab submitting_lab  authors url     title   paper_url       date_submitted
                if( id.endsWith('+') && idx == "taxonomy"){
                    idx = "group"
                }
                if( id.startsWith('EPI') && idx == "taxonomy"){
                    val = `hCoV19/${val}`
                }
                idx = idx.replace("_", " ")
                idx = idx.replace(/^./, idx[0].toUpperCase()); 
                html += `<br><span><strong>${idx}</strong>: ${val}</span>`
            }
        })
        html += "</p></div>"
        return html
    }

    function getRemoteJson(url) {
        return $.ajax({
            type: "GET",
            url: url,
            async: false
        });
    }

    function getIds(obj){
        const g = $(obj).parent().parent()[0]
        const nid = $(g).attr('class').replace('leaf node cid_','')
        const genome_id = assay_stats.responseJSON.tree.nid_to_acc[nid]
        const assay_idx = parseInt($(obj).attr('transform').replace(/[^0-9\-.,]/g, '').split(',')[0])/parseInt($(obj).attr('width'))
        const assay_id = $("text.legend")[assay_idx].innerHTML
        return {'assay_id': assay_id, 'genome_id': genome_id}
    }

    function primerResult(){
        const ids = getIds(this)
        const genome_id = ids.genome_id
        const assay_id = ids.assay_id
        const url = `${assay_result_path}/${assay_id}/${genome_id}.json`
        
        let html = `<div class="mt-2">
        <span><strong>Assay:</strong> ${assay_id}</span><br>
        <span><strong>Target:</strong> ${genome_id}</span>
        </div>`

        const res = getRemoteJson(url)

        // const res = $.getJSON(url)
        if( res.responseText.search("Cannot GET")<0 ){
            $.when(res).done(function(){
                th = res.responseJSON.Thermo
                co = res.responseJSON.Composition
                al = res.responseJSON.Alignments
                va = res.responseJSON.Values
                html += "<div><p>"
                html += `
                <div class="mt-2">
                <span><strong>Amplicon length:</strong> ${co['amplicon length']}bp</span><br>
                <span><strong>Range:</strong> ${co['amplicon range'][0]}..${co['amplicon range'][1]}</span><br>
                </div>
                <div>
                <div class="mt-2">
                <span><strong>Forward-p mismatch:</strong> ${va['Forward Primer']['mismatches']}</span><br>
                <span><strong>Forward-p T<sub>m</sub>:</strong> ${parseFloat(va['Forward Primer']['tm']).toFixed(2)}&deg;C</span>
                </div>
                <div>
                <div class="mt-2">
                <span><strong>Reverse-p mismatch:</strong> ${va['Reverse Primer']['mismatches']}</span><br>
                <span><strong>Reverse-p T<sub>m</sub>:</strong> ${parseFloat(va['Reverse Primer']['tm']).toFixed(2)}&deg;C</span>
                </div>
                <div>
                <div class="mt-2">
                <span><strong>Probe mismatch:</strong> ${va['Probe']['mismatches']}</span><br>
                <span><strong>Probe T<sub>m</sub>:</strong> ${parseFloat(va['Probe']['tm']).toFixed(2)}&deg;C</span><br>
                <span><strong>Probe range:</strong> ${co['probe range'][0]}..${co['probe range'][1]}</span>
                </div>`
                html += "</p><div class='font-italic mb-2'>[Click to display details]<div></div>"
            })  
        }
        else{
            html += "<div class='mt-2'><span>Not available</span><br>"
        }
        return html
    }

    function primerResultDetail(obj){
        const ids = getIds(obj)
        const genome_id = ids.genome_id
        const assay_id = ids.assay_id
        const url = `${assay_result_path}/${assay_id}/${genome_id}.json`
        const res = getRemoteJson(url)

        let title = `Assay: ${assay_id} / Target: ${genome_id}`

        let html = ""

        $.when(res).done(function(){
            th = res.responseJSON.Thermo
            co = res.responseJSON.Composition
            al = res.responseJSON.Alignments
            va = res.responseJSON.Values

            p_pairing = al["Probe"]["pairing"]
            fp_pairing = al["Forward Primer"]["pairing"]
            rp_pairing = al["Reverse Primer"]["pairing"]

            html += `<div>
            <div>
            <span><strong>Amplicon range:</strong> ${co['amplicon range'][0]}..${co['amplicon range'][1]} (${co['amplicon length']}bp)</span><br>
            <span><strong>Target description:</strong> ${res.responseJSON['Common Name']}</span><br>
            <span><strong>Primer GC%:</strong> forward(${parseFloat(co['forward primer %GC ']).toFixed(2)}%), reverse(${parseFloat(co['reverse primer %GC ']).toFixed(2)}%), probe(${parseFloat(co['probe %GC']).toFixed(2)}%)</span><br>
            <span><strong>Min 3' clamp:</strong> ${co["min 3' clamp"]};</span>
            <span><strong>Max 3' clamp:</strong> ${co["max 3' clamp"]}</span><br>
            </div><p>`

            html += `
            <div class="mb-2 font-italic">Assay matching detail</div>
            <table class="table table-striped table-sm">
            <thead>
            <tr>
                <th scope="col">Oligo</th>
                <th scope="col">T<sub>m</sub></th>
                <th scope="col">Hairpin T<sub>m</sub></th>
                <th scope="col">Homodimer T<sub>m</sub></th>
                <th scope="col">Mismatches</th>
                <th scope="col">Gaps</th>
            </tr>
            </thead>
            <tbody class="table-hover">
            <tr>
                <td>Forward</td>
                <td>${parseFloat(va['Forward Primer']['tm']).toFixed(2)}</td>
                <td>${parseFloat(va['Forward Primer']['hairpin tm']).toFixed(2)}</td>
                <td>${parseFloat(va['Forward Primer']['homodimer tm']).toFixed(2)}</td>
                <td>${va['Forward Primer']['mismatches']}</td>
                <td>${va['Forward Primer']['gaps']}</td>
            </tr>
            <tr>
                <td>Reverse</td>
                <td>${parseFloat(va['Reverse Primer']['tm']).toFixed(2)}</td>
                <td>${parseFloat(va['Reverse Primer']['hairpin tm']).toFixed(2)}</td>
                <td>${parseFloat(va['Reverse Primer']['homodimer tm']).toFixed(2)}</td>
                <td>${va['Reverse Primer']['mismatches']}</td>
                <td>${va['Reverse Primer']['gaps']}</td>
            </tr>
            <tr>
                <td>Probe</td>
                <td>${parseFloat(va['Probe']['tm']).toFixed(2)}</td>
                <td>${parseFloat(va['Probe']['hairpin tm']).toFixed(2)}</td>
                <td>${parseFloat(va['Probe']['homodimer tm']).toFixed(2)}</td>
                <td>${va['Probe']['mismatches']}</td>
                <td>${va['Probe']['gaps']}</td>
            </tr>
            </tbody>
            </table>

            <div class="mb-2 font-italic">Alignments</div>
            <table class="table table-striped table-sm table-alignment">
            <thead>
            <tr>
                <th scope="col"> </th>
                <th scope="col">Forward Primer</th>
                <th scope="col">Probe</th>
                <th scope="col">Reverse Primer</th>
                <th scope="col"> </th>
            </tr>
            </thead>
            <tbody class="table-hover">
            <tr>
<td class="text-right">5'-</td>
<td>${al["Forward Primer"]["5'"].replace("5'-","").replace("-3'","")}</td>
<td>${al["Probe"]["5'"].replace("5'-","").replace("-3'","")}</td>
<td>${al["Reverse Primer"]["5'"].replace("5'-","").replace("-3'","")}</td>
<td class="text-left">-3'</td>
            </tr>
            <tr>
<td> </td>
<td>${fp_pairing}</td>
<td>${p_pairing}</td>
<td>${rp_pairing}</td>
<td> </td>
            </tr>
            <tr>
<td class="text-right">3'-</td>
<td>${al["Forward Primer"]["3'"].replace("3'-","").replace("-5'","")}</td>
<td>${al["Probe"]["3'"].replace("3'-","").replace("-5'","")}</td>
<td>${al["Reverse Primer"]["3'"].replace("3'-","").replace("-5'","")}</td>
<td class="text-left">-5'</td>
            </tr>
            </tbody>
            </table>

            <div class="mb-2 font-italic">Thermodynamic information</div>
            <table class="table table-striped table-sm">
            <thead>
            <tr>
                <th scope="col">Oligo</th>
                <th scope="col">dG</th>
                <th scope="col">dH</th>
                <th scope="col">dS</th>
            </tr>
            </thead>
            <tbody class="table-hover">
            <tr>
                <td>Forward primer</td>
                <td>${th['Forward Primer'].dG}</td>
                <td>${th['Forward Primer'].dH}</td>
                <td>${th['Forward Primer'].dS}</td>
            </tr>
            <tr>
                <td>Reverse primer</td>
                <td>${th['Reverse Primer'].dG}</td>
                <td>${th['Reverse Primer'].dH}</td>
                <td>${th['Reverse Primer'].dS}</td>
            </tr>
            <tr>
                <td>Probe</td>
                <td>${th['Probe'].dG}</td>
                <td>${th['Probe'].dH}</td>
                <td>${th['Probe'].dS}</td>
            </tr>
            </tbody>
            </table>`

            html += "</p></div>"
        })

        if(typeof(res.responseJSON)=="undefined"){
            html += "<div>Assay failure</div>"
        }

        return [title, html]
    }

    // Add tooltip
    $('.name').tooltip({
        trigger: 'hover',
        container: 'body',
        delay: 500,
        placement: "left",
        title: nodeMetadata,
        html: true
    });
    $('rect.heatmap').tooltip({
        trigger: 'hover',
        container: 'body',
        delay: 500,
        placement: "bottom",
        title: primerResult,
        boundary: 'scrollParent',
        html: true
    });
    $('rect.heatmap').on('click', function(){
        [title, content] = primerResultDetail(this)
        $('#assayResultModal-title').html(title)
        $('#assayResultModal-body').html(content)
        $('#assayResultModal').modal('show');
    })
}

function loadingMessage(open) {
    if(open){
        $('#spinnerModal').modal('show');
    }
    else{
        $('#spinnerModal').modal('hide');
        $('#spinnerModal').hide()
        $('.loadingSpinner').hide()
    }
}

// Call the dataTables jQuery plugin
$(document).ready(function() {
    loadingMessage(true)

    $('#assayDataTable').DataTable( {
        "pageLength": 5,
        "lengthMenu": [[5, 10, 25, 50, -1], [5, 10, 25, 50, "All"]],
        "ajax": summary_table_json,
        "pagingType": "simple",
        "columns": [
            { "data": "name" },
            { "data": "recall" },
            { "data": "perfect_match" },
            { "data": "1_mm" },
            { "data": "2_mm" },
            { "data": "3_mm_p_fail" }
        ],
        "columnDefs": [
            {   "targets": [0,3,4],
                "class": "nowrap" 
            },
            {
                "targets": 1,
                "render": function (data, type, full) {
                    return (data*100).toString().match(/\d+(\.\d{1,2})?/g)[0];
                }
            },
            {
                "targets": [ 2, 3, 4, 5 ],
                "render": function (data, type, full) {
                    return data.toLocaleString()
                }
            },
            {
                "targets": 0,
                "data": "name",
                "render": function ( data, type, row, meta ) {
                  return `<a class="assay-detail" onclick="assay_mm_charts(this.innerHTML)" style="cursor: pointer;">${data}</a>`;
                }
            }
        ],
    });

    //update spinner
    $.when(assay_stats).done(function(){
        assay_tol_num = Object.keys(assay_stats.responseJSON.tree.assay_stats).length
        collapsed_genome_num = assay_stats.responseJSON.tree.collapsed_genome_num
        $('#spinner-assay-num').html(assay_tol_num.toLocaleString())
        $('#assay-num').html(assay_tol_num.toLocaleString())
        $('#spinner-result-num').html( (assay_tol_num*collapsed_genome_num).toLocaleString() ) 
    })

    $('#assayDataTable_wrapper').removeClass( 'form-inline' );
    $('#sidebarToggle').on("click", function(){
        $('#accordionSidebar').toggleClass("show")
    })

    d3.select("#phyd3").text("Loading phylogeny and heatmap...");
    d3.xml(validation_xml, function(xml) {
        var tree = phyd3.phyloxml.parse(xml);
        phyd3.phylogram.build("#phyd3", tree, phyd3_opts);
        // $('#loading_spinner').toggle();
        $('#phyHeatmap-legend').css("visibility", "visible");
        loadingMessage(false)
    });

    $.when(assay_stats).done(function(){
        updateStats(db_stats.responseJSON, assay_stats.responseJSON)
    })
});