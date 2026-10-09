/* Replace only the two added analysis panels in v015, then save v016. */
(function(){
    var root='/Users/asselism/Documents/Codex/2026-10-03/these-dataset-are-mics-taken-along/outputs/Figure_1_revision';
    var source=root+'/Figure1_v015_revised.ai';
    var log=new File(root+'/tables/illustrator_activity_correction_log.txt');
    function note(s){log.open('a');log.writeln(s);log.close();}
    var prev=app.userInteractionLevel;app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
    try{
        var d=null;
        for(var i=0;i<app.documents.length;i++){
            if(app.documents[i].name==='Figure1_v015_revised.ai'){
                try{if(app.documents[i].fullName.fsName===source){d=app.documents[i];break;}}catch(ignore){}
            }
        }
        if(!d)d=app.open(new File(source));d.activate();
        var layer=d.layers.getByName('Figure 1 additions');
        function replacePanel(stem,x,y,w,h){
            var old=layer.groupItems.getByName(stem);
            var tmp=app.open(new File(root+'/panels/'+stem+'.pdf'));
            var r=tmp.artboards[0].artboardRect;
            tmp.selectObjectsOnActiveArtboard();app.executeMenuCommand('group');
            var g=tmp.selection[0],b=g.geometricBounds;
            var sx=w/(r[2]-r[0]),sy=h/(r[1]-r[3]);
            var copy=g.duplicate(layer,ElementPlacement.PLACEATBEGINNING);
            tmp.close(SaveOptions.DONOTSAVECHANGES);d.activate();
            copy.resize(100*sx,100*sy,true,true,true,true,100*(sx+sy)/2,Transformation.TOPLEFT);
            copy.position=[x+(b[0]-r[0])*sx,-y+(b[1]-r[1])*sy];
            old.remove();copy.name=stem;copy.zOrder(ZOrderMethod.SENDTOBACK);note('Updated '+stem);
        }
        replacePanel('f_inhibition_ridgelines',430,125,185,164);
        replacePanel('g_physicochemical_properties',618,125,257,164);
        d.selection=null;note('Saving v016 with 116 Delta6-only hits.');
        var options=new IllustratorSaveOptions();options.pdfCompatible=true;options.compressed=true;
        d.saveAs(new File(root+'/Figure1_v016_revised.ai'),options);
        note('Saved v016 successfully.');return root+'/Figure1_v016_revised.ai';
    }finally{app.userInteractionLevel=prev;}
})();
