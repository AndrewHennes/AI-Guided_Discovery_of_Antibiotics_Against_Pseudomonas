/* Restore all per-strain hits while preserving the open figure's other edits. */
(function(){
    var root='/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision';
    var out=root+'/panel_g_all_hits_v025';
    var preferred='/Users/asselism/Downloads/Figure1_v024_TPSA_amines_LogP.ai';
    var fallback=root+'/Figure1_v024_TPSA_amines_LogP.ai';
    var target=root+'/Figure1_v025_all_hits_TPSA_amines_LogP.ai';
    if(new File(target).exists)throw new Error('Version 025 already exists; inspect before overwriting.');
    var log=new File(out+'/illustrator_integration_log.txt');
    function note(s){log.open('a');log.writeln(new Date()+' '+s);log.close();}
    function findOpen(path){
        for(var i=0;i<app.documents.length;i++){
            try{if(app.documents[i].fullName.fsName===path)return app.documents[i];}catch(ignore){}
        }
        return null;
    }
    var previous=app.userInteractionLevel;
    app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
    try{
        var d=findOpen(preferred)||findOpen(fallback);
        var wasOpen=(d!==null);
        if(!d)d=app.open(new File(fallback));
        d.activate();note('Source='+d.fullName.fsName+'; was_open='+wasOpen+'; unsaved_changes='+(!d.saved));
        var layer=d.layers.getByName('Figure 1 additions');
        if(layer.placedItems.length!==0)throw new Error('Unexpected placed artwork requires inspection.');
        var old=layer.groupItems.getByName('g_activity_defined_properties');
        var png=new ExportOptionsPNG24();png.antiAliasing=true;png.transparency=false;
        png.artBoardClipping=true;png.horizontalScale=200;png.verticalScale=200;
        d.selection=null;
        d.exportFile(new File(out+'/native_before.png'),ExportType.PNG24,png);
        note('Exported current open artwork before editing.');
        var x=618,top=-125,w=257,h=164;
        var tmp=app.open(new File(out+'/g_TPSA_amines_LogP.pdf'));
        var r=tmp.artboards[0].artboardRect;
        tmp.selectObjectsOnActiveArtboard();app.executeMenuCommand('group');
        var group=tmp.selection[0],bounds=group.geometricBounds;
        var sx=w/(r[2]-r[0]),sy=h/(r[1]-r[3]);
        var copy=group.duplicate(layer,ElementPlacement.PLACEATBEGINNING);
        tmp.close(SaveOptions.DONOTSAVECHANGES);d.activate();
        copy.resize(100*sx,100*sy,true,true,true,true,100*(sx+sy)/2,Transformation.TOPLEFT);
        copy.position=[x+(bounds[0]-r[0])*sx,top+(bounds[1]-r[1])*sy];
        old.remove();copy.name='g_all_hit_properties';copy.zOrder(ZOrderMethod.SENDTOBACK);
        d.selection=null;
        d.exportFile(new File(out+'/native_after.png'),ExportType.PNG24,png);
        note('Replaced only panel g; all-hit counts 8, 5, 39, 124. Saving v025.');
        var options=new IllustratorSaveOptions();options.pdfCompatible=true;options.compressed=true;
        d.saveAs(new File(target),options);
        if(!d.saved)throw new Error('Native save did not complete.');
        note('Saved='+target+'; saved_state='+d.saved+'; placed_items='+layer.placedItems.length);
        if(!wasOpen){d.close(SaveOptions.DONOTSAVECHANGES);note('Closed document opened by this script.');}
        return target;
    }finally{app.userInteractionLevel=previous;}
})();
