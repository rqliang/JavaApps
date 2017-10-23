import org.apache.commons.cli.*;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

public class RunCli {
    public static void main(String[] args) throws Exception {
        Options options = new Options();
        options.addRequiredOption("T", "template_protein",
                true, "template protein in fasta format");
        options.addRequiredOption("F", "DNA_seqs", true,
                "sequencing results in fasta format");
        options.addOption("N", "output filename base", true, "base of alignment result files");
        //options.addOption("P", "out_proteins", true, "output protein sequence file");

        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine cmd = parser.parse(options, args);
            Path templatePath = Paths.get(cmd.getOptionValue("T"));
            Path tpInputFile = Paths.get(cmd.getOptionValue("F"));

            Path pathOutputFile = null;
            Path pathProteinOut = null;
            if(cmd.hasOption("N")) {
                pathOutputFile = Paths.get(cmd.getOptionValue("N") + "_alignment.txt");
                pathProteinOut =  Paths.get(cmd.getOptionValue("N") + "_proteins.fa");
            }





            CliAlignment cliAlignment = new CliAlignment();
            cliAlignment.setTemplateFromFasta(templatePath.toFile());
            cliAlignment.setDnaSequencesFromFasta(tpInputFile.toFile());

            List<ProteinSequence> lp = new ArrayList<ProteinSequence>();
            List<SequencePair<ProteinSequence, AminoAcidCompound>> alignList = new ArrayList<SequencePair<ProteinSequence, AminoAcidCompound>>();
            Map<ProteinSequence, SequencePair<ProteinSequence, AminoAcidCompound>> alignMap = cliAlignment.doAlignmentFromDNASequences();
            for (Map.Entry<ProteinSequence, SequencePair<ProteinSequence, AminoAcidCompound>> entry : alignMap.entrySet()) {
                lp.add(entry.getKey());
                alignList.add(entry.getValue());
            }
            if (pathProteinOut != null) {
                FastaWriterHelper.writeProteinSequence(pathProteinOut.toFile(), (Collection<ProteinSequence>) lp);
            }

            StringBuffer alignPrint = new StringBuffer();
            for (SequencePair<ProteinSequence, AminoAcidCompound> entry : alignList) {
                alignPrint.append("\n----------------------------------------------------\n" + entry.getTarget().getOriginalSequence().getOriginalHeader() + "\n");
                alignPrint.append(entry.toString(60));
            }

            if (pathOutputFile != null) {
                BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(pathOutputFile.toFile()));
                bufferedWriter.write(alignPrint.toString());
                bufferedWriter.close();
            } else {
                System.out.print(alignPrint);
            }

        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("java -jar cli_alignment.jar ", options, true);
            e.printStackTrace();
        }
    }
}