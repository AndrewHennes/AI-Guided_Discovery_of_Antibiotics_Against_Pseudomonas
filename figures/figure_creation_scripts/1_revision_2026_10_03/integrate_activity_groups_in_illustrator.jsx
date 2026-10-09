/* Replace panel g with the three approved property plots and save version 024. */
(function(){
    var root='/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision';
    var source=root+'/Figure1_v023_three_properties.ai';
    var log=new File(root+'/tables/illustrator_activity_group_integration_log.txt');
    function note(s){log.open('a');log.writeln(new Date()+' '+s);log.close();}
    var previous=app.userInteractionLevel;
    app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
    try{
        var d=null;
        for(var i=0;i<app.documents.length;i++){
            if(app.documents[i].name==='Figure1_v023_three_properties.ai'){
                try{if(app.documents[i].fullName.fsName===source){d=app.documents[i];break;}}catch(ignore){}
            }
        }
        if(!d)d=app.open(new File(source));
        d.activate();note('Source '+d.name+'; existing unsaved changes='+(!d.saved));
        var layer=d.layers.getByName('Figure 1 additions');
        if(layer.placedItems.length!==0)throw new Error('Unexpected placed artwork requires inspection.');
        var old=layer.groupItems.getByName('g_three_property_panels');
        // These are the visible page bounds. Clipped marks can extend the group's
        // geometric bounds and must not control placement or scaling.
        var x=618,top=-125,w=257,h=164;
        var tmp=app.open(new File(root+'/panel_g_activity_groups_v024/g_TPSA_amines_LogP.pdf'));
        var r=tmp.artboards[0].artboardRect;
        tmp.selectObjectsOnActiveArtboard();app.executeMenuCommand('group');
        var group=tmp.selection[0],bounds=group.geometricBounds;
        var sx=w/(r[2]-r[0]),sy=h/(r[1]-r[3]);
        var copy=group.duplicate(layer,ElementPlacement.PLACEATBEGINNING);
        tmp.close(SaveOptions.DONOTSAVECHANGES);d.activate();
        copy.resize(100*sx,100*sy,true,true,true,true,100*(sx+sy)/2,Transformation.TOPLEFT);
        copy.position=[x+(bounds[0]-r[0])*sx,top+(bounds[1]-r[1])*sy];
        old.remove();copy.name='g_activity_defined_properties';copy.zOrder(ZOrderMethod.SENDTOBACK);
        // Keep numeric subscripts editable without an unavailable math font.
        for(var ti=0;ti<copy.textFrames.length;ti++){
            var tf=copy.textFrames[ti];
            if(tf.textRange.characterAttributes.textFont.name==='DejaVuSans')
                tf.textRange.characterAttributes.textFont=app.textFonts.getByName('ArialMT');
        }
        d.selection=null;note('Replaced panel g with TPSA, amine count and LogP; group counts 8, 5, 31, 85.');
        var png=new ExportOptionsPNG24();png.antiAliasing=true;png.transparency=false;
        png.artBoardClipping=false;png.horizontalScale=150;png.verticalScale=150;
        d.exportFile(new File(root+'/Figure1_v024_full_artwork.png'),ExportType.PNG24,png);
        note('Exported full artwork preview. Saving v024.');
        var options=new IllustratorSaveOptions();options.pdfCompatible=true;options.compressed=true;
        d.saveAs(new File(root+'/Figure1_v024_TPSA_amines_LogP.ai'),options);
        note('Saved v024. Remaining placed objects='+layer.placedItems.length);
        if(!d.saved)throw new Error('Native save did not complete.');
        // Release this large, fully saved document after the automated edit.
        d.close(SaveOptions.DONOTSAVECHANGES);
        note('Closed the saved document.');
        return root+'/Figure1_v024_TPSA_amines_LogP.ai';
    }finally{app.userInteractionLevel=previous;}
})();
