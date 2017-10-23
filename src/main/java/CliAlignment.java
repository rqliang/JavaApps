import org.apache.commons.cli.*;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;

/**
 * Created by Ruqiang Liang on 4/25/2017.
 */
public class CliAlignment {
    private ProteinSequence template;
    private List<DNASequence> dnaSequences;

    public CliAlignment(){
        template = null;
        dnaSequences = null;
    }
    public CliAlignment(ProteinSequence template, List<DNASequence> dnaSequences) {
        this.template = template;
        this.dnaSequences = dnaSequences;
    }

    public void setTemplateFromFasta(File tp) throws IOException, CompoundNotFoundException {
        //throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        ProteinSequence retP = new ProteinSequence("M");
        LinkedHashMap<String, ProteinSequence> t = FastaReaderHelper.readFastaProteinSequence(tp);
        for (Map.Entry<String, ProteinSequence> entry: t.entrySet() ) {
            String templateName = entry.getValue().getOriginalHeader();
            retP = entry.getValue();
            retP.setOriginalHeader(templateName);
        }
       template = retP;
    }

    public void setDnaSequencesFromFasta(File tp) throws IOException, CompoundNotFoundException {
        DNASequence dnaSequence = new DNASequence("ATGC");
        List<DNASequence> dnaSequenceList = new ArrayList<DNASequence>();
        LinkedHashMap<String, DNASequence> stringDNASequenceLinkedHashMap = FastaReaderHelper.readFastaDNASequence(tp);
        for(Map.Entry<String, DNASequence> entry : stringDNASequenceLinkedHashMap.entrySet()){
            dnaSequenceList.add(entry.getValue());
        }
        dnaSequences = dnaSequenceList;
    }

    private int countStopCodons(ProteinSequence proteinSequence){
        String seq = proteinSequence.toString();
        String s2 = seq.replace("*","");
        return seq.length() - s2.length();
    }
    private List<ProteinSequence> getProteinSequencesFromDNAEntry(DNASequence myDNA) throws CompoundNotFoundException {
        List<ProteinSequence> lp = new ArrayList<ProteinSequence>();
        TranscriptionEngine engine = new TranscriptionEngine.Builder().build();
        Frame[] sixFrames = Frame.getAllFrames();
        Map<Frame, Sequence<AminoAcidCompound>> results = engine.multipleFrameTranslation(myDNA,sixFrames);
        for(Frame frame: sixFrames){
            ProteinSequence proteinSequence = (ProteinSequence)results.get(frame);
            proteinSequence.setOriginalHeader(myDNA.getOriginalHeader() +"_"+ frame.name());
            lp.add(proteinSequence);
        }
        return lp;
    }

    private SequencePair<ProteinSequence, AminoAcidCompound> doAlignmentOneProtein(ProteinSequence prot) {
        SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();
        GapPenalty penalty = new SimpleGapPenalty();
        int gop = 8;
        int extend = 1;
        penalty.setOpenPenalty(gop);
        penalty.setExtensionPenalty(extend);
        PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman;
        smithWaterman = Alignments.getPairwiseAligner(template, prot, Alignments.PairwiseSequenceAlignerType.LOCAL, penalty, matrix);
        SequencePair<ProteinSequence, AminoAcidCompound> pair = smithWaterman.getPair();
        return pair;
    }

    private Map<ProteinSequence, SequencePair<ProteinSequence, AminoAcidCompound>> doAlignmentFromOneDNA(DNASequence dnaSequence) throws CompoundNotFoundException {
        Map<ProteinSequence, SequencePair<ProteinSequence, AminoAcidCompound>> returnMap =
                new HashMap<ProteinSequence, SequencePair<ProteinSequence, AminoAcidCompound>>();
        List<ProteinSequence> proteinSequenceList = getProteinSequencesFromDNAEntry(dnaSequence);
        Map<Integer, SequencePair<ProteinSequence, AminoAcidCompound>> pairMap = new HashMap<Integer, SequencePair<ProteinSequence, AminoAcidCompound>>();
        Integer myMax=0;
        for (ProteinSequence proteinSequence: proteinSequenceList) {
            SequencePair<ProteinSequence, AminoAcidCompound> pair = doAlignmentOneProtein(proteinSequence);
            Integer myIdent = pair.getNumIdenticals();
            pairMap.put(myIdent,pair);
            if (myMax < myIdent) myMax = myIdent;
        }
        SequencePair<ProteinSequence, AminoAcidCompound> thePair = pairMap.get(myMax);
        AlignedSequence<ProteinSequence, AminoAcidCompound> query = thePair.getTarget();
        ProteinSequence theProtein = query.getOriginalSequence();
        theProtein.setDescription(dnaSequence.getOriginalHeader());
        returnMap.put(theProtein, thePair);
        return returnMap;
    }

    public Map<ProteinSequence, SequencePair<ProteinSequence, AminoAcidCompound>> doAlignmentFromDNASequences() throws CompoundNotFoundException {
        Map<ProteinSequence, SequencePair<ProteinSequence, AminoAcidCompound>> returnMap =
                new HashMap<ProteinSequence, SequencePair<ProteinSequence, AminoAcidCompound>>();
        for(DNASequence dnaSequence: dnaSequences) {
            Map<ProteinSequence, SequencePair<ProteinSequence, AminoAcidCompound>> aMap = doAlignmentFromOneDNA(dnaSequence);
            for(Map.Entry<ProteinSequence, SequencePair<ProteinSequence,AminoAcidCompound>> entry: aMap.entrySet()){
                returnMap.put(entry.getKey(),entry.getValue());
            }
        }
        return returnMap;
    }
}
