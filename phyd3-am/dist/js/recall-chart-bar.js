// Set new default font family and font color to mimic Bootstrap's default styling
Chart.defaults.global.defaultFontFamily = 'Arial,sans-serif';
Chart.defaults.global.defaultFontColor = '#858796';

function compare(a, b) {
    // Use toUpperCase() to ignore character casing
    const recallA = a.recall;
    const recallB = b.recall;
    
    let comparison = 0;
    if (recallA > recallB) {
        comparison = -1;
    } else if (recallA < recallB) {
        comparison = 1;
    }
    return comparison;
}

$.when(summary_json).done(function (data) {
    let res = {
        labels: [],
        recalls: []
    }

    let datasets = [{
        label: 'Perfect match',
        backgroundColor: "#506DCE",
        hoverBackgroundColor: "#4e73df",
        data: []
      }, {
        label: '1 mismatch',
        backgroundColor: "#E7CA50",
        data: []
      }, {
        label: '2 mismatch',
        backgroundColor: "#E59244",
        data: []
      }]
  
    //update total genomes and time info
    total_genome_num = data.data[0].perfect_match+data.data[0]['1_mm']+data.data[0]['2_mm_p']+data.data[0].failure
    
    var num_format = (x) =>{
        return (x*100).toString().match(/\d+(\.\d{1,2})?/g)[0]
    }

    $.each(data.data.sort(compare), function (idx, val) {
        if (res.labels.length <= 5) {
            res.labels.push(val.name)
            res.recalls.push(num_format(val.recall))
            res[val.name] = num_format(val.recall)
            datasets[0].data.push(num_format(val.perfect_match/total_genome_num))
            datasets[1].data.push(num_format(val['1_mm']/total_genome_num))
            datasets[2].data.push(num_format(val['2_mm']/total_genome_num))
        }
    })

    // Bar Chart Example
    var ctx = document.getElementById("myBarChart");
    var myBarChart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: res.labels,
            datasets: datasets,
        },
        options: {
            maintainAspectRatio: false,
            layout: {
                padding: {
                    left: 0,
                    right: 10,
                    top: 20,
                    bottom: 0
                }
            },
            hover: {
                onHover: function(e, el) {
                  $("#myBarChart").css("cursor", el[0] ? "pointer" : "default");
                }
            },
            scales: {
            xAxes: [{
                time: {
                unit: '%'
                },
                gridLines: {
                display: false,
                drawBorder: false
                },
                ticks: {
                maxTicksLimit: 6
                },
                maxBarThickness: 25,
                stacked: true,
            }],
            yAxes: [{
                ticks: {
                    min: 98,
                    max: 100,
                    maxTicksLimit: 5,
                    padding: 5,
                    // Include a dollar sign in the ticks
                    //callback: function(value, index, values) {
                    //  return '$' + number_format(value);
                    //}
                },
                stacked: true,
                scaleLabel: {
                    display: true,
                    labelString: "Recall (%)",    
                },
                gridLines: {
                    color: "rgb(234, 236, 244)",
                    zeroLineColor: "rgb(234, 236, 244)",
                    drawBorder: false,
                    borderDash: [2],
                    zeroLineBorderDash: [2]
                }
            }],
            },
            legend: {
                display: false
            },
            tooltips: {
                titleMarginBottom: 10,
                titleFontColor: '#6e707e',
                titleFontSize: 14,
                backgroundColor: "rgb(255,255,255)",
                bodyFontColor: "#858796",
                borderColor: '#dddfeb',
                borderWidth: 1,
                xPadding: 15,
                yPadding: 15,
                caretPadding: 10,
                displayColors: true,
                callbacks: {
                    label: (tooltipItem, data) => {
                        let label = data.datasets[tooltipItem.datasetIndex].label
                        let value = data.datasets[tooltipItem.datasetIndex].data[tooltipItem.index]
                        label = `${label}: ${value}%`
                        return label;
                    },
                    title: (tooltipItem, data) => {
                        let assay_name = tooltipItem[0].xLabel
                        label = `${assay_name} (${res[assay_name]}%)`
                        return label;
                    }
                }
            },
            onClick: function (e) {
                let clickedBar = this.getElementAtEvent(e)[0];
                if (clickedBar) {
                    assay_id = this.data.labels[clickedBar._index];
                    assay_mm_charts(assay_id)
                }
            }
        }
    });
})

function initRes() {
    return {
        labels: [],
        datasets: [{
            label: '1 mismatch',
            backgroundColor: "#858796",
            data: []
          }, {
            label: '2 mismatch',
            backgroundColor: "#76b7b2",
            data: []
          }, {
            label: '3 mismatch',
            backgroundColor: "#edc949",
            data: []
          }, {
            label: '4-10 mismatch',
            backgroundColor: "#f28e2c",
            data: []
          }, {
            label: '11+ mismatch',
            backgroundColor: "#e15759",
            data: []
          }]
    }
}

function assay_mm_charts(assay_id) {
    if(typeof(assay_id)==undefined){
        assay_stats = this.innerHTML
    }
    
    $('#assayStatsModal').modal('show')
    $('#assayStatsModal').on('shown.bs.modal', function(){
        resetMap(assay_id);
        summary_data = {}
        summary_json.responseJSON.data.forEach(element => {
            if(element.name == assay_id){
                summary_data = element
            }
        });
        assay_data = assay_stats.responseJSON.tree.assay_stats[assay_id]

        // update assay description
        fp_seq = assay_stats.responseJSON.tree.assay_stats[assay_id].assay_sequence.forward_primer
        rp_seq = assay_stats.responseJSON.tree.assay_stats[assay_id].assay_sequence.reverse_primer
        pb_seq = assay_stats.responseJSON.tree.assay_stats[assay_id].assay_sequence.probe
        pm_num = summary_data['perfect_match'].toLocaleString()
        desc = `<div class="mt-2 mb-2">
        <span><strong>Assay:</strong> ${assay_id}</span><br>
        <span><strong>Forward primer:</strong> ${fp_seq}</span><br>
        <span><strong>Reverse primer:</strong> ${rp_seq}</span><br>
        <span><strong>Probe:</strong> ${pb_seq}</span><br>
        <span><strong># of perfect match target(s):</strong> ${pm_num}</span><br>
        </div>`

        $('#assayStatsModal-title').html(`${assay_id}: Summary of targets with mismatches`)
        $('#assayResultModal-body div.mm-desc').html(desc)

        // destroy charts before generating new ones
        if(chart1!=false){
            chart1.destroy()
            chart2.destroy()
            chart3.destroy()    
        }

        // chart options
        options = {
            maintainAspectRatio: false,
            responsive: true,
            layout: {
                padding: {
                    left: 10,
                    right: 10,
                    top: 10,
                    bottom: 10
                }
            },
            scales: {
                yAxes: [{
                    scaleLabel: {
                        display: true,
                        labelString: "# of mismatch(es)",    
                    },
                    gridLines: {
                        display: false,
                        drawBorder: true
                    },
                    maxBarThickness: 25,
                    stacked: false,
                }],
                xAxes: [{
                    ticks: {
                        padding: 1,
                        stepSize: 10,
                        maxTicksLimit: 10,
                    },
                    stacked: false,
                    scaleLabel: {
                        display: true,
                        labelString: "# of target(s)",    
                    },
                    gridLines: {
                        color: "rgb(234, 236, 244)",
                        zeroLineColor: "rgb(234, 236, 244)",
                    }
                }],
            },
            legend: {
                display: false
            },
            tooltips: {
                titleMarginBottom: 10,
                titleFontColor: '#6e707e',
                titleFontSize: 14,
                backgroundColor: "rgb(255,255,255)",
                bodyFontColor: "#858796",
                borderColor: '#dddfeb',
                borderWidth: 1,
                xPadding: 15,
                yPadding: 15,
                caretPadding: 10,
                displayColors: true
            }
        }

        // assay-mm-detail
        res = {
            labels: ['1 mismatch', '2 mismatch', '3 mismatch', '4-7 mismatch', '8+ / Failure'],
            datasets: [{
                barPercentage: 0.5,
                barThickness: 6,
                maxBarThickness: 8,
                minBarLength: 2,
                backgroundColor: [
                    "#bab0ab",
                    "#76b7b2",
                    "#edc949",
                    "#f28e2c",
                    "#e15759",
                ], 
                data: [
                    summary_data['1_mm'],
                    summary_data['2_mm'],
                    summary_data['3_mm'],
                    summary_data['4_mm']+summary_data['5_mm']+summary_data['6_mm']+summary_data['7_mm'],
                    parseInt(summary_data['8_mm_p_fail'])
                ]
            }]
        }
    
        ctx = $("#assay-mm-detail")[0]
        ctx.height = 20*res.labels.length+100
        chart1 = new Chart(ctx,{
            type: 'horizontalBar',
            data: {
                labels: res.labels,
                datasets: res.datasets
            },
            options: options
        })

        // data index
        dataset_idx = {
            '1': 0,
            '2': 1,
            '3': 2,
            '4-7': 3,
            '8+': 4
        }

        // time-mm-detail
        options.scales.xAxes[0].stacked = true
        options.scales.yAxes[0].stacked = true
        options.scales.yAxes[0].scaleLabel.labelString = "Time (yyyy-mm)"
        
        res = initRes()
        for (let [month, val] of Object.entries(assay_data.month)) {
            res.labels.push(month)
            for (let [mm_name, idx] of Object.entries(dataset_idx)) {
                res.datasets[idx].data.push(
                    typeof(val[mm_name])=="undefined" ? 0 : val[mm_name]
                )
            }
        }

        ctx = $("#time-mm-detail")[0]
        ctx.height = 20*res.labels.length+100
        chart2 = new Chart(ctx,{
            type: 'horizontalBar',
            data: res,
            options: options        
        })

        // country-mm-detail
        options.scales.yAxes[0].scaleLabel.labelString = "Country"

        function sum(obj) {
            return Object.keys(obj).reduce((sum,key)=>sum+parseFloat(obj[key]||0),0);
        }

        countrySorted = Object.keys(assay_data.country).sort(function(a,b){return sum(assay_data.country[b])-sum(assay_data.country[a])})

        res = initRes()
        // for (let [country, val] of Object.entries(assay_data.country)) {
        countrySorted.forEach((country)=>{
            res.labels.push(country)
            val = assay_data.country[country]
            for (let [mm_name, idx] of Object.entries(dataset_idx)) {
                res.datasets[idx].data.push(
                    typeof(val[mm_name])=="undefined" ? 0 : val[mm_name]
                )
            }
        })

        ctx = $("#country-mm-detail")[0]
        ctx.height = 20*res.labels.length+100
        chart3 = new Chart(ctx,{
            type: 'horizontalBar',
            data: res,
            options: options        
        })
    })
}

let chart1 = false
let chart2 = false
let chart3 = false