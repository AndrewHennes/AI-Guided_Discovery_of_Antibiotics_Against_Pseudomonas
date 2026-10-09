/* Run in Adobe Illustrator. Edits a working copy and saves a new native AI file.
   No command-line arguments. All paths are absolute. */
(function () {
    var root = '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision';
    var source = '/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/work/Figure1_source.ai';
    if (!(new File(source)).exists) {
        var supplied = new File('/Users/asselism/Downloads/Figure1_v014 (3).ai');
        if (!supplied.copy(source)) throw new Error('Unable to create the working copy.');
    }
    var log = new File(root + '/tables/illustrator_build_log.txt');
    function note(s) { log.open('a'); log.writeln(s); log.close(); }
    var previous = app.userInteractionLevel;
    app.userInteractionLevel = UserInteractionLevel.DONTDISPLAYALERTS;
    try {
        var d = null;
        for (var i=0;i<app.documents.length;i++) {
            if (app.documents[i].name === 'Figure1_source.ai') {
                try { if (app.documents[i].fullName.fsName === source) { d=app.documents[i]; break; } } catch (ignored) {}
            }
        }
        if (!d) d=app.open(new File(source));
        d.activate();
        if (d.layers.length > 1 && d.layers[0].name === 'Figure 1 additions') throw new Error('Revision already applied. Open a fresh working copy before rerunning.');
        note('Opened source copy: '+d.name);
        var originalRect=d.artboards[0].artboardRect;
        function selectRegion(rect) {
            d.selection=null;
            var a=d.artboards.add(rect);
            d.artboards.setActiveArtboardIndex(d.artboards.length-1);
            d.selectObjectsOnActiveArtboard();
            var count=d.selection.length;
            a.remove();d.artboards.setActiveArtboardIndex(0);
            return count;
        }
        // Native region selection avoids iterating over hundreds of thousands of markers.
        var count=selectRegion([433,3,650,-123]);
        if (count<10) throw new Error('Count panel region selection was incomplete.');
        app.executeMenuCommand('group');
        var oldCount=d.selection[0];oldCount.name='Original count panel, superseded';
        var archived=d.layers.add();archived.name='Original count panel (hidden)';
        oldCount.move(archived,ElementPlacement.PLACEATBEGINNING);archived.visible=false;
        note('Archived original count panel: '+count+' selected items');
        count=selectRegion([433,-221,875.0155,-469.274]);
        if (count<100) throw new Error('Chemical-space map selection was incomplete.');
        app.executeMenuCommand('group');
        var maps=d.selection[0];maps.name='Original chemical-space maps and legends';maps.translate(0,-70);
        note('Moved map group down 70 pt: '+count+' selected items');
        d.artboards[0].artboardRect=[originalRect[0],originalRect[1],originalRect[2],originalRect[3]-70];
        d.selection=null;
        var added=d.layers.add();added.name='Figure 1 additions';
        function importPDF(stem,x,y,w,h) {
            var temp=app.open(new File(root+'/panels/'+stem+'.pdf'));
            var r=temp.artboards[0].artboardRect;
            temp.selectObjectsOnActiveArtboard();app.executeMenuCommand('group');
            var g=temp.selection[0],b=g.geometricBounds;
            var sx=w/(r[2]-r[0]),sy=h/(r[1]-r[3]);
            var copy=g.duplicate(added,ElementPlacement.PLACEATBEGINNING);
            temp.close(SaveOptions.DONOTSAVECHANGES);d.activate();
            copy.resize(100*sx,100*sy,true,true,true,true,100*(sx+sy)/2,Transformation.TOPLEFT);
            copy.position=[x+(b[0]-r[0])*sx,-y+(b[1]-r[1])*sy];copy.name=stem;
            note('Imported editable PDF: '+stem);
        }
        importPDF('d_updated_compound_counts',430,0,220,123);
        importPDF('f_inhibition_ridgelines',430,125,185,164);
        importPDF('g_physicochemical_properties',618,125,257,164);
        function label(text,x,y) {
            var t=added.textFrames.add();t.contents=text;t.position=[x,-y];
            t.textRange.characterAttributes.textFont=app.textFonts.getByName('Arial-BoldMT');
            t.textRange.characterAttributes.size=13.5;
            var c=new RGBColor();c.red=27;c.green=28;c.blue=24;t.textRange.characterAttributes.fillColor=c;
        }
        label('f',424,125);label('g',607,125);label('h',424,292);label('i',677,292);
        d.selection=null;
        var options=new IllustratorSaveOptions();options.pdfCompatible=true;options.compressed=true;
        note('Saving native AI...');
        d.saveAs(new File(root+'/Figure1_v016_revised.ai'),options);
        note('Saved native AI successfully. Artboard: '+d.artboards[0].artboardRect);
        return root+'/Figure1_v016_revised.ai';
    } finally { app.userInteractionLevel=previous; }
})();
