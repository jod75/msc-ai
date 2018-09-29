// C/C=C(\CCC(C)C1CCC2C1(CCC3C2CCC4=C(C(=O)CC34C)O)C)/C(C)C 

hadoopUrl = 'http://hadoop1:5432'

function sortNumber(a, b) {
    return a - b;
}

function showPage(pageId) {
    document.getElementById('homeContainer').style.display = 'none'
    document.getElementById('toolboxContainer').style.display = 'none'
    document.getElementById('resultsContainer').style.display = 'none'

    document.getElementById(pageId).style.display = 'block'
}

function getMoleculeDataFromSmiles(smiles) {
    require(['https://www.lactame.com/lib/openchemlib/5.2.0/openchemlib-full.js'],
        function(OCL) {
            try {
                var molecule = OCL.Molecule.fromSmiles(smiles);
                var mf = molecule.getMolecularFormula();
                var properties = new OCL.MoleculeProperties(molecule);
                var result = {
                    //svg: molecule.toSVG(400, 350),
                    mw: mf.relativeWeight,
                    em: mf.absoluteWeight,
                    mf: mf.formula,
                    logP: properties.logP,
                    logS: properties.logS,
                    psa: properties.polarSurfaceArea,
                    error: ''
                };
            } catch (err) {
                result = {
                    //svg: '',
                    mw: 0,
                    em: 0,
                    mf: '',
                    logP: 0,
                    logS: 0,
                    psa: 0,
                    error: err
                }
            }

            // update data on screen
            showMoelculeInformation(result);

            // unload module in case of error as subsequent calls return same error
            if (result.error != "") {
                require.undef('https://www.lactame.com/lib/openchemlib/5.2.0/openchemlib-full.js');
            }
        });

    // render query SVG
    $.ajax({
        url: hadoopUrl + "getSmilesSVG/" + encodeURIComponent(smiles) + "/smiles",
        type: "get",
        datatype: "json",
        success: function(response) {
            document.getElementById('smilesCanvas').innerHTML = response;
        },
        error: function() {
            var img = document.createElement("IMG");
            img.src = "no_image.png";
            document.getElementById('smilesCanvas').appendChild(img);
        }
    });

    document.getElementById("moleculeRegNo").innerText = "Checking ChEMBL data...";
    var img = document.createElement("IMG");
    img.src = "ajax-loader.gif";
    document.getElementById("moleculeRegNo").appendChild(img);
    $.ajax({
        url: hadoopUrl + "isLigandInChEMBL/" + encodeURIComponent(smiles),
        type: "get",
        datatype: "json",
        success: function(response) {
            response = response.replace(/([\[\]])/g, "")
            if (response != "") {
                document.getElementById("moleculeRegNo").innerText = "ChEMBL (molRegNo): " + response;
                document.getElementById('moleculeBindings').innerText = "Checking if there are any recorded target proteins in local ChEMBL database.";
                var img = document.createElement("IMG");
                img.src = "ajax-loader.gif";
                document.getElementById("moleculeBindings").appendChild(img);
                $.ajax({
                    url: hadoopUrl + "getLigandBindings/" + encodeURIComponent(smiles),
                    type: "get",
                    datatype: "json",
                    success: function(response) {
                        var targets = JSON.parse(response);
                        if (targets.length == 0) {
                            document.getElementById('moleculeBindings').innerText = "No known protein targets in local ChEMBL database.";
                        } else {
                            document.getElementById('moleculeBindings').innerText = "Binds to proteins: " + targets[0][1];
                        }
                        for (var i = 1; i < targets.length; i++) {
                            document.getElementById('moleculeBindings').innerText += (", " + targets[i][1]);
                        }
                    }
                });
            } else {
                document.getElementById("moleculeRegNo").innerText = "Not available in local ChEMBL database";
                document.getElementById('moleculeBindings').innerText = "";
            }
        },
        error: function() {
            showConnectionAlert();
            document.getElementById("moleculeRegNo").innerText = "Cannot connect to ICS5200 server.";
        }
    });
}

function showMoelculeInformation(info) {
    var moleculeDataField = document.getElementById('moleculeData');
    if (info.error == '') {
        document.getElementById('virtualScreen').style.display = 'block';
        moleculeDataField.innerHTML = ("Relative weight: " + info.mw + " g/mol<br>");
        moleculeDataField.innerHTML += ("Absolute weight: " + info.em + " g/mol<br>");
        moleculeDataField.innerHTML += ("Formula: " + info.mf.replace(/(\d+)/g, '<sub>$1</sub>') + "<br>");
        moleculeDataField.innerHTML += ("logP: " + info.logP + "<br>");
        moleculeDataField.innerHTML += ("logS: " + info.logS + "<br>");
        moleculeDataField.innerHTML += ("Polar surface area: " + info.psa + ' &#8491;<sup>2</sup>');
    } else {
        //document.getElementById('toolboxChemInfo').style.display = 'none';
        moleculeDataField.innerHTML = info.error;
    }
}

function drawMolecule() {
    var smiles = document.getElementById('txtToolboxSearchLigand').value.trim();
    getMoleculeDataFromSmiles(smiles);
}

function selectProteinToolbox() {
    document.getElementById('btnToolboxSearchProtein').removeAttribute('disabled', false);
    document.getElementById('txtToolboxSearchProtein').removeAttribute('disabled', false);
    document.getElementById('btnToolboxSearchLigand').setAttribute('disabled', true);
    document.getElementById('txtToolboxSearchLigand').setAttribute('disabled', true);
    document.getElementById('cardToolboxLigand').style.display = 'none';
    document.getElementById('cardToolboxProtein').style.display = 'block';
}

function selectLigandToolbox() {
    document.getElementById('btnToolboxSearchProtein').setAttribute('disabled', true);
    document.getElementById('txtToolboxSearchProtein').setAttribute('disabled', true);
    document.getElementById('btnToolboxSearchLigand').removeAttribute('disabled', false);
    document.getElementById('txtToolboxSearchLigand').removeAttribute('disabled', false);
    document.getElementById('cardToolboxLigand').style.display = 'block';
    document.getElementById('cardToolboxProtein').style.display = 'none';
}

function showConnectionAlert() {
    $("#noconnection-alert").fadeTo(20000, 500).slideUp(500, function() {
        $("#noconnection-alert").slideUp(500);
    });
}

$(document).ready(function() {
    // ---- Hide alerts
    $("#noconnection-alert").hide();

    // ---- Set containers
    showPage('toolboxContainer');

    // ---- By default enable ligand search
    document.getElementById('rdToolboxSearchLigand').checked = true;
    document.getElementById('cardToolboxProtein').style.display = 'none';
    drawMolecule();

    // ---- Format molecule similarity table
    var molSimTable = $('#ligandSimilarityTable').DataTable({
        colReorder: true,
        dom: 'Bfrtip',
        buttons: [
            'copyHtml5',
            'excelHtml5',
            'csvHtml5',
            'pdfHtml5'
        ]
    });
    molSimTable.column(2).visible(false);
    molSimTable.column(3).visible(false);
    molSimTable.column(4).visible(false);
    molSimTable.column(10).visible(false);
    molSimTable.column(13).visible(false);
    molSimTable.column(14).visible(false);
    molSimTable.column(11).visible(false);
    molSimTable.colReorder.order([0, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 11, 13, 14, 15]);

    // ------------------------------------------------------------------------------------------------------------
    // Nav menu

    // ---- Callback functions

    $('#toolboxMenuItem').click(function() {
        showPage('toolboxContainer');
    });

    $('#homeMenuItem').click(function() {
        showPage('homeContainer');
    });

    // ------------------------------------------------------------------------------------------------------------
    // Toolbox

    $('#rdToolboxSearchProtein').click(function() {
        selectProteinToolbox();
    });

    $('#rdToolboxSearchLigand').click(function() {
        selectLigandToolbox();
    });

    $('#btnToolboxSearchLigand').click(drawMolecule);

    $('#btnToolboxSearchProtein').click(function() {
        var sequence = document.getElementById("txtToolboxSearchProtein").value.trim();

        $.ajax({
            url: hadoopUrl + "isProteinInChEMBL/" + sequence,
            type: "get",
            datatype: "json",
            success: function(response) {
                var proteins = JSON.parse(response);
                if (proteins.length > 0) {
                    document.getElementById("toolboxProteinDescription").innerText = "Protein detail: ";
                    for (var i = 0; i < proteins.length; i++) {
                        document.getElementById("toolboxProteinDescription").innerHTML += proteins[i][1] + ", " + proteins[i][0] + "<br>";
                    }
                } else {
                    document.getElementById("toolboxProteinDescription").innerText = "Protein not found in local ChEMBL database.";
                }
            },
            error: function() {
                showConnectionAlert();
            }
        });
    });

    $('#btnToolboxScreenLigands').click(function() {
        var smiles = encodeURIComponent(document.getElementById('txtToolboxSearchLigand').value);
        var fp = "Morgan";
        if (document.getElementById('ckToolboxMACCSnFingerprint').checked) {
            fp = "MACCS";
        }
        var sim = "Tanimoto";
        var th = document.getElementById("txtToolboxTanimotoCoefficient").value;

        document.getElementById("queryLigandMolRegNo").innerText = document.getElementById('moleculeRegNo').innerText;

        // render query SVG
        $.ajax({
            url: hadoopUrl + "getSmilesSVG/" + smiles + "/smiles",
            type: "get",
            datatype: "json",
            success: function(response) {
                document.getElementById('queryLigandSVG').innerHTML = response;
            },
            error: function() {
                document.getElementById('queryLigandSVG').src = "no_image.png"
            }
        });

        // run experiment using smiles
        $.ajax({
            url: hadoopUrl + "doLigandExperimentFromSmiles/" + smiles + "/" + fp + "/" + sim + "/" + th,
            type: "get",
            datatype: "json",
            success: function(response) {
                var res = JSON.parse(response);
                molSimTable.rows().clear();
                molSimTable.colReorder.reset();
                molSimTable.rows.add(res).draw();
                molSimTable.colReorder.order([0, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 11, 13, 14, 15]);

                var knownLigandsUniqueCompId = document.getElementById("knownLigandsUniqueCompId");
                knownLigandsUniqueCompId.innerHTML = "( ";
                var uniqueCompIds = molSimTable.column(14).data().unique().toArray().sort(sortNumber);
                for (var i = 0; i < uniqueCompIds.length; i++) {
                    knownLigandsUniqueCompId.innerHTML += (uniqueCompIds[i] + " ");
                }
                knownLigandsUniqueCompId.innerHTML += ")";
                showPage('resultsContainer');
            }
        });
    });

    $('#btnToolboxScreenResultsLigands').click(function() {
        showPage('resultsContainer');
    });


    // ------------------------------------------------------------------------------------------------------------
    // Ligand results

    $('#ligandSimilarityTable tbody').on('click', 'tr', function() {
        if (molSimTable.row(this).length > 0) {
            document.getElementById('knownLigandMolRegNo').innerText = molSimTable.row(this).data()[0];
            $.ajax({
                url: hadoopUrl + "getSmilesSVG/" + molSimTable.row(this).data()[0] + "/mol",
                type: "get",
                datatype: "json",
                success: function(response) {
                    document.getElementById('knownLigandSVG').innerHTML = response
                }
            });

            $.ajax({
                url: hadoopUrl + "getSmiles/" + molSimTable.row(this).data()[0],
                type: "get",
                datatype: "json",
                success: function(response) {
                    document.getElementById('knownLigandSmiles').innerHTML = response.replace(/"/g, "")
                }
            });
        }
    });

    // ------------------------------------------------------------------------------------------------------------
    // Ajax spinner

    // http://gasparesganga.com/labs/jquery-loading-overlay/
    $(document).ajaxStart(function() {
        $.LoadingOverlay("show");
    });

    $(document).ajaxStop(function() {
        $.LoadingOverlay("hide");
    });
});