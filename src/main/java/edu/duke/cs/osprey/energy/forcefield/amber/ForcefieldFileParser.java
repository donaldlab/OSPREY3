/*
 ** This file is part of OSPREY 3.0
 **
 ** OSPREY Protein Redesign Software Version 3.0
 ** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
 **
 ** OSPREY is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License version 2
 ** as published by the Free Software Foundation.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
 **
 ** OSPREY relies on grants for its development, and since visibility
 ** in the scientific literature is essential for our success, we
 ** ask that users of OSPREY cite our papers. See the CITING_OSPREY
 ** document in this distribution for more information.
 **
 ** Contact Info:
 **    Bruce Donald
 **    Duke University
 **    Department of Computer Science
 **    Levine Science Research Center (LSRC)
 **    Durham
 **    NC 27708-0129
 **    USA
 **    e-mail: www.cs.duke.edu/brd/
 **
 ** <signature of Bruce Donald>, Mar 1, 2018
 ** Bruce Donald, Professor of Computer Science
 */

package edu.duke.cs.osprey.energy.forcefield.amber;

import com.beust.jcommander.internal.Lists;
import one.util.streamex.StreamEx;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.nio.file.Files;
import java.nio.file.Path;
import java.io.IOException;
import java.util.Optional;
import java.util.stream.Collectors;

public class ForcefieldFileParser {

    private final InputStream parmFile;

    private record TwoTuple<T>(T head, T tail) {
    }

    // For use in dihedrals where If IPT .eq. 'X ' .and. LPT .eq. 'X ' then
    // any dihedrals in the system involving the atoms "JPT" and "KPT" are assigned the same parameters.
    public static final AtomSymbolAndMass WildcardAtom = new AtomSymbolAndMass("X", 0);
    public static final AtomSymbolAndMass UnmatchedAtom = new AtomSymbolAndMass("?", 0);
    private final Path frcmod;

    public ForcefieldFileParser(InputStream parmFile) {
        this(parmFile, null);
    }

    public ForcefieldFileParser(InputStream parmFile, Path frcmod) {
        this.parmFile = parmFile;
        this.frcmod = frcmod;
    }

    private String title;
    private final List<AtomSymbolAndMass> atomSymbolsAndMasses = Lists.newArrayList();
    private final List<AtomSymbolAndMass> hydrophilicAtoms = Lists.newArrayList();
    private final List<BondLengthParameter> bondLengthParameters = Lists.newArrayList();
    private final List<BondAngleParameter> bondAngleParameters = Lists.newArrayList();
    private final List<DihederalParameter> dihederalParameters = Lists.newArrayList();
    private final List<ImproperDihederalParameter> improperDihederalParameters = Lists.newArrayList();
    private final List<HBond10_12PotentialParameter> hbond10_12PotentialParameters = Lists.newArrayList();
    private final List<EquivalencingAtom> equivalencingAtomsForNonBonded6_12PotentialParameters = Lists.newArrayList();
    private final List<SlaterKirkwoodParameter> slaterKirkwoodParameters = Lists.newArrayList();
    private final List<VanDerWaalsRadius> vanDerWaalsRadii = Lists.newArrayList();
    private final List<Six12PotentialCoefficient> six12PotentialCoefficients = Lists.newArrayList();

    public List<AtomSymbolAndMass> atomSymbolsAndMasses() {
        return new ArrayList<>(atomSymbolsAndMasses);
    }

    private static <T extends HasAtoms> List<T> filterOutUnmatchedAtom(List<T> lst) {
        return StreamEx.of(lst).filter(x -> !x.atoms().contains(ForcefieldFileParser.UnmatchedAtom)).toList();
    }

    public List<AtomSymbolAndMass> hydrophilicAtoms() {
        return StreamEx.of(hydrophilicAtoms).filter(x -> x != UnmatchedAtom).toList();
    }

    public List<BondLengthParameter> bondLengthParameters() {
        return filterOutUnmatchedAtom(bondLengthParameters);
    }

    public List<BondAngleParameter> bondAngleParameters() {
        return filterOutUnmatchedAtom(bondAngleParameters);
    }

    public List<DihederalParameter> dihederalParameters() {
        return filterOutUnmatchedAtom(dihederalParameters);
    }

    public List<ImproperDihederalParameter> improperDihederalParameters() {
        return filterOutUnmatchedAtom(improperDihederalParameters);
    }

    public List<HBond10_12PotentialParameter> hbond10_12PotentialParameters() {
        return filterOutUnmatchedAtom(hbond10_12PotentialParameters);
    }

    public List<EquivalencingAtom> equivalencingAtomsForNonBonded6_12PotentialParameters() {
        return filterOutUnmatchedAtom(equivalencingAtomsForNonBonded6_12PotentialParameters);
    }

    public List<SlaterKirkwoodParameter> slaterKirkwoodParameters() {
        return filterOutUnmatchedAtom(slaterKirkwoodParameters);
    }

    public List<VanDerWaalsRadius> vanDerWaalsRadii() {
        return filterOutUnmatchedAtom(vanDerWaalsRadii);
    }

    public List<Six12PotentialCoefficient> six12PotentialCoefficients() {
        return filterOutUnmatchedAtom(six12PotentialCoefficients);
    }


    /**
     * This method reads all of the parameters from an Amber parm.dat file into the fields. If a frcmod file is
     * provided, it reads those, too.
     * @throws IOException Thrown if there's an issue reading the parm.dat or frcmod file
     */
    public void read() throws IOException {
        var parmLines = new BufferedReader(new InputStreamReader(parmFile, StandardCharsets.UTF_8))
                .lines()
                .collect(Collectors.toList());

        title = parmLines.get(0);
        parmLines = parmLines.subList(1, parmLines.size());

        var secondBlock = takeUntilBlankLine(parmLines);
        var atomSymbolsAndMassesLines = secondBlock.head;
        var frcmodAtomSymbols = getFromFrcmod("MASS");
        atomSymbolsAndMasses.addAll(
            concatWithReplacement(
                    parseAtomSymbolsAndMasses(atomSymbolsAndMassesLines),
                    parseAtomSymbolsAndMasses(frcmodAtomSymbols.orElse(List.of())),
                    (a, b) -> a.KNDSYM().equals(b.KNDSYM())
            )
        );

        var thirdBlock = secondBlock.tail;
        var hydrophilicAtomLines = thirdBlock.get(0);
        hydrophilicAtoms.addAll(parseHydrophilicAtoms(hydrophilicAtomLines));

        var fourthBlock = takeUntilBlankLine(secondBlock.tail.subList(1, secondBlock.tail.size()));
        var bondLengthLines = fourthBlock.head;
        var frcmodBondLengthLines = getFromFrcmod("BOND");
        bondLengthParameters.addAll(
                concatWithReplacement(
                        parseBondLengthParameters(bondLengthLines),
                        parseBondLengthParameters(frcmodBondLengthLines.orElse(List.of())),
                        (a, b) -> pairwiseEquals(
                                List.of(a.IBT(), a.JBT()),
                                List.of(b.IBT(), b.JBT())
                        )
                )
        );

        var fifthBlock = takeUntilBlankLine(fourthBlock.tail);
        var bondAngleLines = fifthBlock.head;
        var frcmodBondAngleLines = getFromFrcmod("ANGL");
        bondAngleParameters.addAll(
                concatWithReplacement(
                        parseBondAngleParameters(bondAngleLines),
                        parseBondAngleParameters(frcmodBondAngleLines.orElse(List.of())),
                        (a, b) -> pairwiseEquals(
                                List.of(a.ITT(), a.JTT(), a.KTT()),
                                List.of(b.ITT(), b.JTT(), b.KTT()))

                )
        );

        var sixthBlock = takeUntilBlankLine(fifthBlock.tail);
        var dihedralLines = sixthBlock.head;
        var frcmodDihedralLines = getFromFrcmod("DIHE");
        dihederalParameters.addAll(
                concatWithReplacement(
                        parseDihedralParameters(dihedralLines),
                        parseDihedralParameters(frcmodDihedralLines.orElse(List.of())),
                        (a, b) -> pairwiseEquals(
                                List.of(a.IPT(), a.JPT(), a.KPT(), a.LPT()),
                                List.of(b.IPT(), b.JPT(), b.KPT(), b.LPT()))
                )
        );

        var seventhBlock = takeUntilBlankLine(sixthBlock.tail);
        var improperDihedralLines = seventhBlock.head;
        var frcmodImproperDihedralLines = getFromFrcmod("IMPR");
        improperDihederalParameters.addAll(
                concatWithReplacement(
                        parseImproperDihedralParameters(improperDihedralLines),
                        parseImproperDihedralParameters(frcmodImproperDihedralLines.orElse(List.of())),
                        (a, b) -> pairwiseEquals(
                                List.of(a.IPT(), a.JPT(), a.KPT(), a.LPT()),
                                List.of(b.IPT(), b.JPT(), b.KPT(), b.LPT())
                        )
                )
        );

        var eighthBlock = takeUntilBlankLine(seventhBlock.tail);
        var hBond10_12PotentialLines = eighthBlock.head;
        var frcmodHBond10_12Lines = getFromFrcmod("HBON");
        hbond10_12PotentialParameters.addAll(
                concatWithReplacement(
                        parseHBond10_12PotentialParameters(hBond10_12PotentialLines)        ,
                        parseHBond10_12PotentialParameters(frcmodHBond10_12Lines.orElse(List.of())),
                        (a, b) -> pairwiseEquals(
                                List.of(a.KT1(), a.KT2()),
                                List.of(b.KT1(), b.KT2())
                        )
                )
        );

        var ninthBlock = takeUntilBlankLine(eighthBlock.tail);
        var equivalencingAtomsForNonBonded6_12PotentialParametersLines = ninthBlock.head;
        equivalencingAtomsForNonBonded6_12PotentialParameters.addAll(parseEquivalencingAtomfForNonBonded6_12PotentialParameters(equivalencingAtomsForNonBonded6_12PotentialParametersLines));

        var tenthBlock = ninthBlock.tail;
        var labelKindLine = tenthBlock.get(0);
        var six12ParameterLines = takeUntilBlankLine(tenthBlock.subList(1, tenthBlock.size()));
        var frcmod6_12ParameterLines = getFromFrcmod("NONB").orElse(List.of());
        String six12ParameterType = labelKindLine.substring(10, 12);

        switch (six12ParameterType) {
            case "SK" -> slaterKirkwoodParameters.addAll(
                    concatWithReplacement(
                            parseSlaterKirkwoodParameters(six12ParameterLines.head),
                            parseSlaterKirkwoodParameters(frcmod6_12ParameterLines),
                            (a, b) -> a.LTYNB().equals(b.LTYNB())
                    )
            );
            case "RE" -> vanDerWaalsRadii.addAll(
                    concatWithReplacement(
                            parseVanDerWaalsRadiiParameters(six12ParameterLines.head),
                            parseVanDerWaalsRadiiParameters(frcmod6_12ParameterLines),
                            (a, b) -> a.LTYNB().equals(b.LTYNB())
                    )
            );
            case "AC" -> six12PotentialCoefficients.addAll(
                    concatWithReplacement(
                            parse6_12PotentialCoefficientParameters(six12ParameterLines.head),
                            parse6_12PotentialCoefficientParameters(frcmod6_12ParameterLines),
                            (a, b) -> a.LTYNB().equals(b.LTYNB())
                    )
            );
            default -> throw new RuntimeException(String.format("The 6-12 potential parameter kind %s is not known", six12ParameterType));
        }
    }

	private static String eolSubstring(String orig, int start, int end) {
		var len = Math.min(orig.length(), end);
		return orig.substring(start, len);
	}

    private Optional<List<String>> getFromFrcmod(String kind) throws IOException {
        if (frcmod == null) {
            return Optional.empty();
        }

        // the head line is a remark line
        var frcmodLines = StreamEx.of(Files.readAllLines(frcmod))
                .skip(1)
                .toList();

        if (frcmodLines.isEmpty()) {
            return Optional.empty();
        }

        return findBlock(takeUntilBlankLine(frcmodLines), kind);
    }

    private Optional<List<String>> findBlock(TwoTuple<List<String>> block, String targetKind) {
        var foundKind = block.head.get(0).substring(0, 4);

        if (foundKind.equals(targetKind)) {
            return Optional.of(block.head.subList(1, block.head.size()));
        }

        if (block.tail.stream().anyMatch(l -> !l.isBlank())) {
            block = takeUntilBlankLine(block.tail);
            return findBlock(block, targetKind);
        }

        return Optional.empty();
    }

    private interface AreSame<T> {
        boolean check(T a, T b);
    }

    private <T> List<T> concatWithReplacement(List<T> head, List<T> tail, AreSame<T> determiner) {
        return StreamEx.of(
                head.stream()
                     .filter(headEl -> tail
                        .stream()
                        .noneMatch(tailEl -> determiner.check(headEl, tailEl))))
                .append(tail.stream()).toList();
    }

    private <T> boolean pairwiseEquals(List<T> a, List<T> b) {
        return StreamEx.zip(a, b, Object::equals).allMatch(j -> j);
    }

    private Six12PotentialCoefficient parse6_12PotentialCoefficientParameter(String line) {
        // FORMAT(2X,A2,6X,3F10.6)
        return new Six12PotentialCoefficient(
                findAtomSymbolFromName(line.substring(2, 4)),
                Float.parseFloat(line.substring(10, 20)),
                Float.parseFloat(line.substring(20, 30))
        );
    }

    private List<Six12PotentialCoefficient> parse6_12PotentialCoefficientParameters(List<String> lines) {
        return StreamEx.of(lines)
                .map(this::parse6_12PotentialCoefficientParameter)
                .toList();
    }

    private VanDerWaalsRadius parseVanDerWaalsRadiusParameter(String line) {
        // FORMAT(2X,A2,6X,3F10.6)
        return new VanDerWaalsRadius(
                findAtomSymbolFromName(line.substring(2, 4)),
                Float.parseFloat(line.substring(10, 20)),
                Float.parseFloat(eolSubstring(line, 20, 30))
        );
    }

    private List<VanDerWaalsRadius> parseVanDerWaalsRadiiParameters(List<String> lines) {
        return StreamEx.of(lines)
                .map(this::parseVanDerWaalsRadiusParameter)
                .toList();
    }

    private SlaterKirkwoodParameter parseSlaterKirkwoodParamter(String line) {
        // FORMAT(2X,A2,6X,3F10.6)
        return new SlaterKirkwoodParameter(
                findAtomSymbolFromName(line.substring(2, 4)),
                Float.parseFloat(line.substring(10, 20)),
                Float.parseFloat(line.substring(20, 30)),
                Float.parseFloat(line.substring(30, 40))
        );
    }

    private List<SlaterKirkwoodParameter> parseSlaterKirkwoodParameters(List<String> lines) {
        return StreamEx.of(lines)
                .map(this::parseSlaterKirkwoodParamter)
                .toList();
    }

    private List<EquivalencingAtom> parseEquivalencingAtomfForNonBonded6_12PotentialParameters(List<String> lines) {
        return StreamEx.of(lines)
                .flatMap(l -> Arrays.stream(
                        l.substring(4).split(" {2}"))
                        .map(rest -> new EquivalencingAtom(
                                findAtomSymbolFromName(l.substring(0, 2)), findAtomSymbolFromName(rest))
                        ))
                .filter(x -> x.IEQV() != UnmatchedAtom && x.IORG() != UnmatchedAtom)
                .distinct(EquivalencingAtom::IEQV)
                .toList();
    }

    private HBond10_12PotentialParameter parseHBond10_12PotentialParameter(String line) {
        return new HBond10_12PotentialParameter(
                findAtomSymbolFromName(line.substring(2, 4)),
                findAtomSymbolFromName(line.substring(6, 8)),
                Float.parseFloat(line.substring(8, 18)),
                Float.parseFloat(line.substring(18, 28))
        );
    }

    private List<HBond10_12PotentialParameter> parseHBond10_12PotentialParameters(List<String> lines) {
        return StreamEx.of(lines)
                .map(this::parseHBond10_12PotentialParameter)
                .toList();
    }

    private ImproperDihederalParameter parseImproperDihedralParameter(String line) {
        return new ImproperDihederalParameter(
                findAtomSymbolOrWildcardFromName(line.substring(0, 2)),
                findAtomSymbolOrWildcardFromName(line.substring(3, 5)),
                findAtomSymbolOrWildcardFromName(line.substring(6, 8)),
                findAtomSymbolOrWildcardFromName(line.substring(9, 11)),
                Float.parseFloat(line.substring(11, 26)),
                Float.parseFloat(line.substring(26, 41)),
                Float.parseFloat(line.substring(41, Math.min(56, line.length()))) // sometimes the last float isn't given the standard 15 chars
        );
    }

    private List<ImproperDihederalParameter> parseImproperDihedralParameters(List<String> lines) {
        return StreamEx.of(lines)
                .map(this::parseImproperDihedralParameter)
                .toList();
    }

    private DihederalParameter parseDihedralParameter(String line) {
        // FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
        return new DihederalParameter(
                findAtomSymbolOrWildcardFromName(line.substring(0, 2)),
                findAtomSymbolOrWildcardFromName(line.substring(3, 5)),
                findAtomSymbolOrWildcardFromName(line.substring(6, 8)),
                findAtomSymbolOrWildcardFromName(line.substring(9, 11)),
                Integer.parseInt(line.substring(11, 15).trim()),
                Float.parseFloat(line.substring(15, 30)),
                Float.parseFloat(line.substring(30, 45)),
                Float.parseFloat(line.substring(45, Math.min(60, line.length()))) // sometimes the last float isn't given the standard 15 chars
        );
    }

    private List<DihederalParameter> parseDihedralParameters(List<String> lines) {
        return StreamEx.of(lines)
                .map(this::parseDihedralParameter)
                .toList();
    }

    private BondAngleParameter parseBondAngleParameter(String line) {
        // FORMAT(A2,1X,A2,1X,A2,2F10.2)
        return new BondAngleParameter(
                findAtomSymbolFromName(line.substring(0, 2)),
                findAtomSymbolFromName(line.substring(3, 5)),
                findAtomSymbolFromName(line.substring(6, 8)),
                Float.parseFloat(line.substring(8, 18)),
                Float.parseFloat(line.substring(18, 28)));
    }

    private List<BondAngleParameter> parseBondAngleParameters(List<String> lines) {
        return StreamEx.of(lines)
                .map(this::parseBondAngleParameter)
                .toList();
    }

    private BondLengthParameter parseBondLengthParameter(String line) {
        return new BondLengthParameter(
                findAtomSymbolFromName(line.substring(0, 2)),
                findAtomSymbolFromName(line.substring(3, 5)),
                Float.parseFloat(line.substring(5, 15)),
                Float.parseFloat(eolSubstring(line, 15, 25)));
    }

    private List<BondLengthParameter> parseBondLengthParameters(List<String> lines) {
        return StreamEx.of(lines)
                .map(this::parseBondLengthParameter)
                .toList();
    }

    private AtomSymbolAndMass findAtomSymbolOrWildcardFromName(String s) {
        return s.trim().equals("X") ? WildcardAtom : findAtomSymbolFromName(s);
    }

    private AtomSymbolAndMass findAtomSymbolFromName(String s) {
        return StreamEx.of(atomSymbolsAndMasses)
                .filter(a -> a.KNDSYM().equals(s.trim()))
                .findFirst()
                .orElse(UnmatchedAtom);
    }

    private static TwoTuple<List<String>> takeUntilBlankLine(List<String> lines) {
        var beforeBreak = StreamEx.of(lines)
                .takeWhile(l -> !l.trim().isEmpty())
                .toList();

        var afterBreak = StreamEx.of(lines)
                .skip(beforeBreak.size() + 1)
                .toList();

        return new TwoTuple<>(beforeBreak, afterBreak);
    }

    private List<AtomSymbolAndMass> parseHydrophilicAtoms(String line) {
        // FORMAT(20(A2,2X))
        return StreamEx.of(line.split(" {2}"))
                .map(this::findAtomSymbolFromName)
                .filter(x -> x != UnmatchedAtom)
                .toList();
    }

    private static AtomSymbolAndMass parseAtomSymbolAndMass(String line) {
        // FORMAT(A2,2X,F10.2x,f10.2)
        // in actuality, the lines appear to be FORMAT(A2,1X,F10.2x,f10.2)
        return new AtomSymbolAndMass(line.substring(0, 2).trim(), Float.parseFloat(line.substring(3, 13)));
    }

    private List<AtomSymbolAndMass> parseAtomSymbolsAndMasses(List<String> lines) {
        return StreamEx.of(lines)
                .map(ForcefieldFileParser::parseAtomSymbolAndMass)
                .toList();
    }
}
