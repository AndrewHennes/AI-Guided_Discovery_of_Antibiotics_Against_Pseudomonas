/* Remove only the two inspected stray placed figures and save version 017. */
(function(){
    var root='/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision';
    var d;try{d=app.documents.getByName('Figure1_v016_revised.ai');}
    catch(missingDocument){d=app.open(new File(root+'/Figure1_v016_revised.ai'));}
    d.activate();
    var layer=d.layers.getByName('Figure 1 additions');
    var ab=d.artboards[0].artboardRect;
    var targets=[];
    for(var i=0;i<layer.placedItems.length;i++){
        var p=layer.placedItems[i],b=p.visibleBounds;
        var above=b[1]>ab[1]+100 && b[3]>=ab[1]-1;
        var below=b[3]<ab[3]-200 && b[1]<ab[3]+100;
        if(above||below)targets.push(p);
    }
    if(targets.length!==2)throw new Error('Expected the two inspected stray placed figures; found '+targets.length+'.');
    var report=new File(root+'/tables/removed_stray_artwork.txt');report.open('w');
    report.writeln('Original artboard\t'+ab);
    for(var j=targets.length-1;j>=0;j--){
        var q=targets[j],source='Unavailable';try{source=q.file.fsName;}catch(ignore){}
        report.writeln('Removed placed artwork\t'+q.visibleBounds+'\t'+source);
        q.remove();
    }
    report.writeln('Remaining placed items in additions layer\t'+layer.placedItems.length);
    report.writeln('Cleaned document bounds\t'+d.visibleBounds);report.close();
    d.selection=null;
    var png=new ExportOptionsPNG24();png.antiAliasing=true;png.transparency=false;png.artBoardClipping=false;png.horizontalScale=150;png.verticalScale=150;
    d.exportFile(new File(root+'/Figure1_v017_full_artwork.png'),ExportType.PNG24,png);
    var options=new IllustratorSaveOptions();options.pdfCompatible=true;options.compressed=true;
    d.saveAs(new File(root+'/Figure1_v017_cleaned.ai'),options);
    return root+'/Figure1_v017_cleaned.ai';
})();
