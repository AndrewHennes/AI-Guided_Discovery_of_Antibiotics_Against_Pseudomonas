(function () {
    var directory = new File($.fileName).parent;
    var stem = directory.fsName + "/full_figure_2/deepflux_predictions_and_hit_properties_figure_2";
    var input = new File(stem + ".pdf");
    if (!input.exists) {
        throw new Error("Run reproduce_figure_2.py before saving the Illustrator document.");
    }
    var document = app.open(input);
    var options = new IllustratorSaveOptions();
    options.pdfCompatible = true;
    options.compressed = true;
    options.compatibility = Compatibility.ILLUSTRATOR17;
    document.saveAs(new File(stem + ".ai"), options);
})();
