//init
var map_assay_id = undefined;
var mismatch_type = "8";
var mismatch_text = 'Failure';
var collection_start_date = moment().subtract(29, 'days');
var collection_end_date = moment();

var marker_scale = 10;
var marker_colors = {
    "0": "#4e73df",
    "1": "#858796",
    "2": "#76b7b2",
    "3": "#edc949",
    "5": "#f28e2c",
    "8": "#e15759",
    "A": "#db4655",
}
var marker_group;

//mismatch type select2
var mismatch_types = [
    {
        id: 0,
        text: 'Perfect Match'
    },
    {
        id: 1,
        text: '1 Mismatch'
    },
    {
        id: 2,
        text: '2 Mismatches'
    },
    {
        id: 3,
        text: '3 Mismatches'
    },
    {
        id: 5,
        text: '4-7 Mismatches'
    },
    {
        id: 8,
        text: '8+/Failure'
    },
    {
        id: 'A',
        text: 'Total failures'
    }
];

var mapMismatchSelect2 = $("#map-mismatch-select2").select2({
    data: mismatch_types,
    minimumResultsForSearch: -1,
    multiple: false,
    selectOnClose: true,
    width: '100%'
});

$('#map-mismatch-select2').on("change", function (e) {
    $('#map-mismatch-select2 option:selected').each(function () {
        mismatch_type = $(this).val();
        mismatch_text = $(this).text();
    });
    //update map markers
    setMapMarkers();
});

//target collection date range
function setDates(start, end) {
    $('#collection_date_range span').html(start.format('MMMM D, YYYY') + ' - ' + end.format('MMMM D, YYYY'));
    collection_start_date = start;
    collection_end_date = end;
    //update map markers
    setMapMarkers();
}

$('#collection_date_range').daterangepicker({
    startDate: collection_start_date,
    endDate: collection_end_date,
    ranges: {
        'Today': [moment(), moment()],
        'Yesterday': [moment().subtract(1, 'days'), moment().subtract(1, 'days')],
        'Last 7 Days': [moment().subtract(6, 'days'), moment()],
        'Last 30 Days': [moment().subtract(29, 'days'), moment()],
        'Last 90 Days': [moment().subtract(89, 'days'), moment()],
        'Last 180 Days': [moment().subtract(179, 'days'), moment()],
    }
}, setDates);

//leaflet map
var map_mismatch = L.map('map-mismatch').setView([51.505, -0.09], 13);

function resetMap(assay_id) {
    map_assay_id = undefined;
    //reset mismatch select
    mapMismatchSelect2.val('A').trigger("change");
    //reset target collection date range
    collection_start_date = moment().subtract(180, 'days');
    collection_end_date = moment();
    var drp = $('#collection_date_range').data('daterangepicker');
    drp.setStartDate(collection_start_date);
    drp.setEndDate(collection_end_date);
    setDates(collection_start_date, collection_end_date);

    //map
    //leaflet map
    if (map_mismatch) {
        map_mismatch.off();
        map_mismatch.remove();
    }
    map_mismatch = L.map('map-mismatch').setView([35, 23], 2);
    // L.tileLayer('https://tiles.stadiamaps.com/tiles/outdoors/{z}/{x}/{y}{r}.png', {
    //     maxZoom: 20,
    //     minZoom: 2,
    //     attribution: '&copy; <a href="https://stadiamaps.com/">Stadia Maps</a>, &copy; <a href="https://openmaptiles.org/">OpenMapTiles</a> &copy; <a href="http://openstreetmap.org">OpenStreetMap</a> contributors'
    // }).addTo(map_mismatch);

    // L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    //     maxZoom: 20,
    //     minZoom: 2,
    //     attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
    // }).addTo(map_mismatch);
    
    L.tileLayer('http://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png', {
        maxZoom: 20,
        minZoom: 2,
        attribution: ''
    }).addTo(map_mismatch);

    // L.tileLayer('https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png', {
    //     maxZoom: 20,
    //     minZoom: 2,
    //     attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>'
    // }).addTo(map_mismatch);

    /*Legend specific*/
    var legend = L.control({ position: "topright" });
    legend.onAdd = function (map) {
        var div = L.DomUtil.create("div", "map-legend");
        div.innerHTML += '<i style="background: ' + marker_colors['0'] + '"></i><span>Perfect Match</span><br>';
        div.innerHTML += '<i style="background: ' + marker_colors['1'] + '"></i><span>1 Mismatch</span><br>';
        div.innerHTML += '<i style="background: ' + marker_colors['2'] + '"></i><span>2 Mismatches</span><br>';
        div.innerHTML += '<i style="background: ' + marker_colors['3'] + '"></i><span>3 Mismatches</span><br>';
        div.innerHTML += '<i style="background: ' + marker_colors['5'] + '"></i><span>4-7 Mismatches</span><br>';
        div.innerHTML += '<i style="background: ' + marker_colors['8'] + '"></i><span>8+/Failures</span><br>';
        div.innerHTML += '<i style="background: ' + marker_colors['A'] + '"></i><span>Total failures(3+)</span><br>';

        return div;
    };
    legend.addTo(map_mismatch);

    //markers
    map_assay_id = assay_id;
    setMapMarkers();

    //re-adjust the width/height bounds of the L.Map's container
    map_mismatch.invalidateSize();
}

function setMapMarkers() {
    if (!map_assay_id) {
        return;
    }
    //get data based on the mismatch type and collection date range
    var map_data = getMapData();
    var latlngs = country_latlngs.responseJSON;

    /* delete old markers */
    if (marker_group) {
        map_mismatch.removeLayer(marker_group);
    }
    marker_group = L.layerGroup();

    //add new markers
    map_data.map(item => {
        var marker = L.circleMarker([latlngs[item.country]['latitude'], latlngs[item.country]['longitude']], {
            color: marker_colors[mismatch_type],
            fillColor: marker_colors[mismatch_type],
            fillOpacity: 0.5,
            radius: Math.log(item.count + 1) * marker_scale
        });
        marker.bindPopup(map_assay_id + "<br>Country: " + latlngs[item.country]['country'] + "<br>" + mismatch_text + ": " + item.count);
        marker.on('mouseover', function (e) {
            this.openPopup();
        });
        marker.on('mouseout', function (e) {
            this.closePopup();
        });
        marker_group.addLayer(marker);
    });
    map_mismatch.addLayer(marker_group);
}

function getMapData() {
    var data = assay_val_results.responseJSON[map_assay_id];
    var start_date = moment(collection_start_date).format('YYYY-MM-DD');
    var end_date = moment(collection_end_date).format('YYYY-MM-DD');

    var map_data = {};
    if (mismatch_type === '0') {
        map_data = getMapMarkerData(data['Perfect match'], start_date, end_date);
    }
    else if (mismatch_type === '1') {
        map_data = getMapMarkerData(data['1 mismatch'], start_date, end_date);
    }
    else if (mismatch_type === '2') {
        map_data = getMapMarkerData(data['2 mismatches'], start_date, end_date);
    }
    else if (mismatch_type === '3') {
        map_data = getMapMarkerData(data['3 mismatches'], start_date, end_date);
    }
    else if (mismatch_type === '5') {
        map_data = getMapMarkerData(data['4-7 mismatches'], start_date, end_date);
    }
    else if (mismatch_type === '8') {
        map_data = getMapMarkerData(data['8+/failures'], start_date, end_date);
    }
    else if (mismatch_type === 'A') {
        map_data = getMapMarkerData(data['Total failures'], start_date, end_date);
    }

    return map_data;
}

//get a list of countries and sum the number of targets with collection date within the selected date range
function getMapMarkerData(data, start_date, end_date) {
    var marker_countries = [];
    var collection_date;
    if (data) {
        Object.keys(data).map(country => {
            var targets = data[country];
            var total = 0;
            Object.keys(targets).map(key => {
                collection_date = moment(key).format('YYYY-MM-DD');
                if (key >= start_date && key <= end_date) {
                    total += targets[key];
                }
            });
            if (total > 0) {
                marker_countries.push({ country: country, count: total });
            }
        });
    }

    return marker_countries;
}