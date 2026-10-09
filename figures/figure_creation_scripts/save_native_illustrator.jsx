(function () {
    var directory = new File($.fileName).parent;
    var stem = directory.fsName + "/full_figure/antipseudomonal_screening_and_chemical_space_figure_1";
    var input = new File(stem + ".pdf");
    if (!input.exists) {
        throw new Error("Run reproduce_figure_1.py before saving the Illustrator document.");
    }
    var document = app.open(input);
    var options = new IllustratorSaveOptions();
    options.pdfCompatible = true;
    options.compressed = true;
    options.compatibility = Compatibility.ILLUSTRATOR17;
    document.saveAs(new File(stem + ".ai"), options);
})();
