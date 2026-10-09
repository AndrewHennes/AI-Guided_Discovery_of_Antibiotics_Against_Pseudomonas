/* Preserve current edits in the open cleaned document, replace f and g, save v018. */
(function(){
    var root='/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision';
    var source=root+'/Figure1_v017_cleaned.ai';
    var log=new File(root+'/tables/illustrator_violin_update_log.txt');
    function note(s){log.open('a');log.writeln(new Date()+' '+s);log.close();}
    var previous=app.userInteractionLevel;
    app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
    try{
        var d=null;
        for(var i=0;i<app.documents.length;i++){
            if(app.documents[i].name==='Figure1_v017_cleaned.ai'){
                try{if(app.documents[i].fullName.fsName===source){d=app.documents[i];break;}}catch(ignore){}
            }
        }
        if(!d)d=app.open(new File(source));
        d.activate();
        note('Source '+d.name+'; existing unsaved changes='+(!d.saved));
        var layer=d.layers.getByName('Figure 1 additions');
        if(layer.placedItems.length!==0)throw new Error('Unexpected placed artwork requires inspection.');
        // Confirm both targets exist before replacing either one.
        layer.groupItems.getByName('f_inhibition_ridgelines');
        layer.groupItems.getByName('g_physicochemical_properties');
        function replacePanel(oldName,newName){
            var old=layer.groupItems.getByName(oldName),target=old.geometricBounds;
            var x=target[0],top=target[1],w=target[2]-target[0],h=target[1]-target[3];
            var tmp=app.open(new File(root+'/panels/'+newName+'.pdf'));
            var r=tmp.artboards[0].artboardRect;
            tmp.selectObjectsOnActiveArtboard();app.executeMenuCommand('group');
            var group=tmp.selection[0],bounds=group.geometricBounds;
            var sx=w/(r[2]-r[0]),sy=h/(r[1]-r[3]);
            var copy=group.duplicate(layer,ElementPlacement.PLACEATBEGINNING);
            tmp.close(SaveOptions.DONOTSAVECHANGES);d.activate();
            copy.resize(100*sx,100*sy,true,true,true,true,100*(sx+sy)/2,Transformation.TOPLEFT);
            copy.position=[x+(bounds[0]-r[0])*sx,top+(bounds[1]-r[1])*sy];
            old.remove();copy.name=newName;copy.zOrder(ZOrderMethod.SENDTOBACK);
            note('Replaced '+oldName+' with '+newName+' at '+target);
        }
        replacePanel('f_inhibition_ridgelines','f_inhibition_log_density_violins');
        replacePanel('g_physicochemical_properties','g_physicochemical_properties');
        d.selection=null;
        var png=new ExportOptionsPNG24();png.antiAliasing=true;png.transparency=false;
        png.artBoardClipping=false;png.horizontalScale=150;png.verticalScale=150;
        d.exportFile(new File(root+'/Figure1_v018_full_artwork.png'),ExportType.PNG24,png);
        note('Exported full artwork preview. Saving v018.');
        var options=new IllustratorSaveOptions();options.pdfCompatible=true;options.compressed=true;
        d.saveAs(new File(root+'/Figure1_v018_violins.ai'),options);
        note('Saved v018. Remaining placed objects='+layer.placedItems.length);
        return root+'/Figure1_v018_violins.ai';
    }finally{app.userInteractionLevel=previous;}
})();
