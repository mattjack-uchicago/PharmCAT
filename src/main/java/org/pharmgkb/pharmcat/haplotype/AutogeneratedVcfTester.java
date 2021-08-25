package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Sets;
import org.apache.commons.lang3.time.DurationFormatUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.IoUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * This class runs all autogenerated test VCFs against the {@link NamedAlleleMatcher}.
 * <ul>
 *   <li>By default, test will pass if expected result is one of the highest matching results from the
 *     {@link NamedAlleleMatcher}.</li>
 *   <li>In fuzzy-match mode, test will pass if expected result matches any result from the
 *     {@link NamedAlleleMatcher}.</li>
 *   <li>In exact-match mode, test will only pass if the {@link NamedAlleleMatcher} produces a single result that
 *     matches the expected result.</li>
 * </ul>
 *
 * @author Mark Woon
 */
public class AutogeneratedVcfTester implements AutoCloseable {
  private static final ResultSerializer sf_resultSerializer = new ResultSerializer();
  private final Path m_exemptionsFile;
  private final Path m_outputDir;
  private final boolean m_saveData;
  private final boolean m_exactMatchOnly;
  private final boolean m_fuzzyMatch;
  private int m_numTests;
  private int m_numWarnings;
  private int m_numFailures;
  private final ErrorWriter m_warningWriter;
  private final ErrorWriter m_failureWriter;


  private AutogeneratedVcfTester(Path outputDir, boolean saveData, boolean exactMatchOnly, boolean fuzzyMatch)
      throws IOException {
    Preconditions.checkArgument(!(exactMatchOnly && fuzzyMatch));
    m_exemptionsFile = DataManager.DEFAULT_DEFINITION_DIR.resolve(DataManager.EXEMPTIONS_JSON_FILE_NAME);
    if (!Files.isRegularFile(m_exemptionsFile)) {
      throw new IllegalStateException("Cannot find exemptions file: " + m_exemptionsFile);
    }
    m_outputDir = outputDir;
    m_warningWriter = new ErrorWriter(m_outputDir.resolve("autogenerated_test_report_warnings.txt"));
    m_failureWriter = new ErrorWriter(m_outputDir.resolve("autogenerated_test_report.txt"));
    m_saveData = saveData;
    m_exactMatchOnly = exactMatchOnly;
    m_fuzzyMatch = fuzzyMatch;
  }

  @Override
  public void close() {
    IoUtils.closeQuietly(m_warningWriter);
    IoUtils.closeQuietly(m_failureWriter);
  }


  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("vcf", "vcf-dir", "test VCF directory", true, "vcf")
          .addOption("o", "output-dir", "output directory", true, "o")
          .addOption("g", "gene", "restrict to specific gene", false, "g")
          .addOption("s", "save", "save all results")
          .addOption("e", "exact-match-only", "only pass if matcher produces single exact match")
          .addOption("f", "fuzzy-match", "pass if matcher produces any match")
          ;

      cliHelper.execute(args, cli -> {
        try{
          Path vcfDir = cliHelper.getValidDirectory("vcf", false);
          Path outputDir = cliHelper.getValidDirectory("o", true);
          boolean exact = cliHelper.hasOption("e");
          boolean fuzzy = cliHelper.hasOption("f");
          if (exact && fuzzy) {
            System.out.println("exact-match-only and fuzzy-match are mutually exclusive");
            return 1;
          }

          try (AutogeneratedVcfTester tester =
                   new AutogeneratedVcfTester(outputDir, cliHelper.hasOption("s"), exact, fuzzy)) {
            Stopwatch stopwatch = Stopwatch.createStarted();

            if (cli.hasOption("g")) {
              Path geneDir = vcfDir.resolve(Objects.requireNonNull(cliHelper.getValue("g")));
              if (Files.isDirectory(geneDir)) {
                System.out.println(geneDir + " is not a valid directory");
              }
              tester.testGene(geneDir);
            } else {
              tester.testAllGenes(vcfDir);
            }

            stopwatch.stop();
            System.out.println("Done.");
            System.out.println(DurationFormatUtils.formatDurationHMS(stopwatch.elapsed(TimeUnit.MILLISECONDS)));
          }
          return 0;

        } catch (Exception ex) {
          ex.printStackTrace();
          return 1;
        }
      });
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  private void testAllGenes(Path vcfDir) throws IOException {
    Stopwatch stopwatch = Stopwatch.createStarted();

    try (DirectoryStream<Path> geneDirStream = Files.newDirectoryStream(vcfDir)) {
      for (Path geneDir : geneDirStream) {
        // TODO: remove after v1 is released
        // not performing matching on CYP2D6 for now
        if (geneDir.getFileName().toString().equals("CYP2D6")) {
          continue;
        }
        if (Files.isDirectory(geneDir)) {
          testGene(geneDir);
        }
      }
    }

    String elapsedTime = DurationFormatUtils.formatDurationHMS(stopwatch.elapsed(TimeUnit.MILLISECONDS));
    m_warningWriter.printSummary(m_numTests, m_numWarnings, m_numFailures, elapsedTime);
    m_failureWriter.printSummary(m_numTests, m_numWarnings, m_numFailures, elapsedTime);
  }

  private void testGene(Path geneDir) throws IOException {

    String gene = geneDir.getFileName().toString();
    System.out.println("Testing " + gene + "...\n");
    Path definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve(gene + "_translation.json");
    if (!Files.isRegularFile(definitionFile)) {
      throw new IllegalStateException("Cannot find definition file for " + gene + ": " + definitionFile);
    }
    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(definitionFile);
    definitionReader.readExemptions(m_exemptionsFile);
    DefinitionExemption exemptions = definitionReader.getExemption(gene);
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(definitionReader, true, false);

    m_warningWriter.startGene(gene);
    m_failureWriter.startGene(gene);
    try (DirectoryStream<Path> vcfStream = Files.newDirectoryStream(geneDir)) {
      for (Path vcfFile : vcfStream) {
        if (Files.isRegularFile(vcfFile) && vcfFile.toString().endsWith(".vcf")) {
          try {
            test(gene, namedAlleleMatcher, vcfFile, definitionReader, exemptions);
          } catch (RuntimeException ex) {
            throw new RuntimeException("Error on " + vcfFile, ex);
          }
        }
      }
    }
    System.out.println();
  }


  private void test(String gene, NamedAlleleMatcher namedAlleleMatcher, Path vcfFile, DefinitionReader definitionReader,
      @Nullable DefinitionExemption exemption) throws IOException {

    m_numTests += 1;
    VcfReader vcfReader = new VcfReader(vcfFile);
    String expectedDiplotype = vcfReader.getVcfMetadata().getRawProperties().get("PharmCATnamedAlleles").get(0);
    List<String> expectedAlleles = Arrays.asList(expectedDiplotype.split("/"));
    boolean hasUnknownCall = expectedAlleles.contains("?");
    boolean hasComboCall = !hasUnknownCall && vcfFile.getFileName().toString().contains("noCall");

    Result result = namedAlleleMatcher.call(vcfFile);
    Set<DiplotypeMatch> matches = result.getGeneCalls().get(0).getDiplotypes();
    boolean gotExpected = false;
    // only compare top candidates
    List<DiplotypeMatch> topPairs = new ArrayList<>();
    List<DiplotypeMatch> alternatePairs = new ArrayList<>();
    if (matches.size() > 0) {
      int topScore = matches.iterator().next().getScore();
      for (DiplotypeMatch match : matches) {
        if (match.getName().equals(expectedDiplotype)) {
          gotExpected = true;
        }
        if (match.getScore() == topScore) {
          topPairs.add(match);
        } else {
          alternatePairs.add(match);
        }
      }
    }

    if (hasUnknownCall || hasComboCall) {
      if (topPairs.size() > 0) {
        fail(vcfFile, result, topPairs, alternatePairs, "no call (" + expectedDiplotype + ")", null, exemption);
      }
      return;
    }
    if (m_exactMatchOnly && matches.size() > 1) {
      String extraMsg = null;
      if (gotExpected) {
        List<String> errors = checkOverlaps(definitionReader.getPositions(gene), matches);
        if (errors.size() > 0) {
          if (errors.size() == 1) {
            extraMsg = errors.get(0);
          } else {
            extraMsg = String.join("\n  ", errors);
          }
        } else {
          warn(vcfFile, result, topPairs, alternatePairs, expectedDiplotype, null, exemption);
          return;
        }
      }
      fail(vcfFile, result, topPairs, alternatePairs, expectedDiplotype, extraMsg, exemption);
      return;
    }
    if (topPairs.size() != 1 && !m_fuzzyMatch) {
      fail(vcfFile, result, topPairs, alternatePairs, expectedDiplotype, null, exemption);
      return;
    }

    Collections.sort(expectedAlleles);
    String expected = String.join("/", expectedAlleles);

    if (!m_fuzzyMatch) {
      if (!isMatch(expected, topPairs.get(0))) {
        fail(vcfFile, result, topPairs, alternatePairs, expectedDiplotype, null, exemption);
      }

    } else {
      boolean foundMatch = result.getGeneCalls().get(0).getDiplotypes().stream()
          .anyMatch(dm -> isMatch(expected, dm));
      if (!foundMatch) {
        fail(vcfFile, result, topPairs, alternatePairs, expectedDiplotype, null, exemption);
      }
    }
    if (m_saveData) {
      saveData(vcfFile, result);
    }
  }

  private List<String> checkOverlaps(VariantLocus[] positions, Set<DiplotypeMatch> matches) {
    //noinspection UnstableApiUsage
    Set<Set<DiplotypeMatch>> combinations = Sets.combinations(matches, 2);
    List<String> overlaps = new ArrayList<>();
    for (Set<DiplotypeMatch> combo : combinations) {
      Iterator<DiplotypeMatch> it = combo.iterator();
      DiplotypeMatch m1 = it.next();
      DiplotypeMatch m2 = it.next();
      int size = m1.getHaplotype1().getHaplotype().getAlleles().length;
      String[] alleles1 = new String[size];
      String[] alleles2 = new String[size];
      for (int x = 0; x < size; x += 1) {
        alleles1[x] = buildAllele(m1, x);
        alleles2[x] = buildAllele(m2, x);
      }
      if (!Arrays.equals(alleles1, alleles2)) {
        StringBuilder errBuilder = new StringBuilder();
        for (int x = 0; x < size; x += 1) {
          if (!Objects.requireNonNull(alleles1[x]).equals(alleles2[x])) {
            if (errBuilder.length() > 0) {
              errBuilder.append("\n");
            }
            errBuilder.append("  Mismatch in ")
                .append(positions[x])
                .append(": ")
                .append(alleles1[x])
                .append(" vs. ")
                .append(alleles2[x]);
          }
        }
        overlaps.add(m1.getName() + " and " + m2.getName() + " DO NOT overlap:\n" + errBuilder);
      }
    }
    return overlaps;
  }

  private String buildAllele(DiplotypeMatch m, int x) {
    SortedSet<String> bases = new TreeSet<>();
    if (m.getHaplotype1().getHaplotype().getAlleles()[x] != null) {
      String allele = m.getHaplotype1().getHaplotype().getAlleles()[x];
      if (allele.length() == 1) {
        bases.addAll(Iupac.lookup(allele).getBases());
      } else {
        bases.add(allele);
      }
    }
    if (m.getHaplotype2().getHaplotype().getAlleles()[x] != null) {
      String allele = m.getHaplotype2().getHaplotype().getAlleles()[x];
      if (allele.length() == 1) {
        bases.addAll(Iupac.lookup(allele).getBases());
      } else {
        bases.add(allele);
      }
    }
    if (bases.size() > 0) {
      return String.join(",", bases);
    }
    return null;
  }


  private boolean isMatch(String expected, DiplotypeMatch match) {
    List<String> rezAlleles = Arrays.asList(match.getName().split("/"));
    Collections.sort(rezAlleles);
    String rez = String.join("/", rezAlleles);
    return rez.equals(expected);
  }


  private void warn(Path vcfFile, Result result, List<DiplotypeMatch> topPairs, List<DiplotypeMatch> alternatePairs,
      String expected, String extraMsg, @Nullable DefinitionExemption exemption) throws IOException {
    addError(vcfFile, result, topPairs, alternatePairs, expected, extraMsg, exemption, true);
  }

  private void fail(Path vcfFile, Result result, List<DiplotypeMatch> topPairs, List<DiplotypeMatch> alternatePairs,
      String expected, String extraMsg, @Nullable DefinitionExemption exemption) throws IOException {
    addError(vcfFile, result, topPairs, alternatePairs, expected, extraMsg, exemption, false);
  }

  private void addError(Path vcfFile, Result result, List<DiplotypeMatch> topPairs, List<DiplotypeMatch> alternatePairs,
      String expected, String extraMsg, @Nullable DefinitionExemption exemption, boolean warn) throws IOException {
    String type = " [FAILURE]";
    if (warn) {
      type = " [WARNING]";
      m_numWarnings += 1;
    } else {
      m_numFailures += 1;
    }
    //noinspection UnstableApiUsage
    String baseFilename = com.google.common.io.Files.getNameWithoutExtension(vcfFile.getFileName().toString());
    String actual = topPairs.stream()
        .map(DiplotypeMatch::getName)
        .collect(Collectors.joining(", "));
    String alt = null;
    if (alternatePairs.size() > 0) {
      alt = alternatePairs.stream()
          .map(m -> m.getName() + " (" + m.getScore() + ")")
          .collect(Collectors.joining(", "));
    }

    System.out.println("* " + baseFilename + type);
    System.out.println("  Expected: " + expected);
    System.out.println("    Actual: " + actual);
    if (alt != null) {
      System.out.println("      Alts: " + alt);
    }
    if (extraMsg != null) {
      System.out.println("  " + extraMsg);
    }
    System.out.println();

    StringBuilder errBuilder = new StringBuilder()
        .append(baseFilename).append(type)
        .append("\n")
        .append("  Expected: ").append(expected)
        .append("\n")
        .append("    Actual: ").append(actual);
    if (topPairs.size() > 0) {
      errBuilder.append(" (")
          .append(topPairs.get(0).getScore())
          .append(")");
    }
    errBuilder.append("\n");
    if (alt != null) {
      errBuilder.append("      Alts: ")
          .append(alt)
          .append("\n");
    }
    if (extraMsg != null) {
      errBuilder.append(extraMsg)
          .append("\n");
    }
    if (exemption != null) {
      if (exemption.isAssumeReference() == Boolean.FALSE) {
        errBuilder.append("EXEMPTION: not assuming reference")
            .append("\n");
      }
      if (exemption.isAllHits() == Boolean.FALSE) {
        errBuilder.append("EXEMPTION: top candidates only")
            .append("\n");
      }
    }
    errBuilder.append("\n");
    if (warn) {
      m_warningWriter.println(errBuilder.toString());
    } else {
      m_failureWriter.println(errBuilder.toString());
    }

    saveData(vcfFile, result);
  }

  private void saveData(Path vcfFile, Result result) throws IOException {
    //noinspection UnstableApiUsage
    String baseFilename = com.google.common.io.Files.getNameWithoutExtension(vcfFile.getFileName().toString());
    Files.copy(vcfFile, m_outputDir.resolve(vcfFile.getFileName()), StandardCopyOption.REPLACE_EXISTING);
    sf_resultSerializer.toJson(result, m_outputDir.resolve(baseFilename + ".json"));
    sf_resultSerializer.toHtml(result, m_outputDir.resolve(baseFilename + ".html"));
  }


  private static class ErrorWriter implements AutoCloseable {
    private final PrintWriter m_writer;
    private String m_gene;
    private boolean m_started;
    private boolean m_geneStarted;


    private ErrorWriter(Path file) throws IOException {
      m_writer = new PrintWriter(Files.newBufferedWriter(file));
    }


    @Override
    public void close() {
      IoUtils.closeQuietly(m_writer);
    }

    private void startGene(String gene) {
      m_gene = gene;
      m_geneStarted = false;
      m_writer.flush();
    }

    private void println(String msg) {
      if (!m_geneStarted) {
        m_geneStarted = true;
        if (m_started) {
          m_writer.println();
          m_writer.println("------------------------------");
          m_writer.println();
        } else {
          m_started = true;
        }
        m_writer.println(m_gene);
        m_writer.println();
      }
      m_writer.println(msg);
    }

    private void printSummary(int numTests, int numWarnings, int numFailures, String elapsedTime) {
      NumberFormat numFormatter = NumberFormat.getInstance();
      m_writer.println("");
      m_writer.println("# tests    = " + numFormatter.format(numTests));
      m_writer.println("# passed   = " + numFormatter.format(numTests - numFailures - numWarnings));
      m_writer.println("# warnings = " + numFormatter.format(numWarnings));
      m_writer.println("# failed   = " + numFormatter.format(numFailures));
      m_writer.println("Elapsed    = " + elapsedTime);
    }
  }
}
