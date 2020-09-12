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

package edu.duke.cs.osprey.markstar.visualizer;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;
import javafx.collections.ObservableList;
import javafx.geometry.Point2D;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.paint.Color;
import javafx.scene.shape.*;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.scene.transform.Transform;

import java.io.*;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.IntStream;

import static edu.duke.cs.osprey.tools.Log.log;

public class KStarTreeNode implements Comparable<KStarTreeNode>{
    public static final Pattern p = Pattern.compile("((~\\+)*)\\(([^)]+)\\)->\\((.*)\\)\\: ?\\[(.*)\\]->\\[(.*)\\](.*)");
    protected static final boolean debug = false;
    protected static Random[] colorSeeds;
    protected BigDecimal overallUpperBound;
    public int level = -1;
    protected String[] assignments;
    protected int[] confAssignments;
    protected BigDecimal upperBound;
    protected BigDecimal lowerBound;
    List<KStarTreeNode> children = new ArrayList<>();
    protected KStarTreeNode parent;
    protected boolean visible =  false;
    protected boolean expanded = false;
    protected boolean childrenRendered = false;
    protected static Text statText;
    protected Group bandGroup;
    protected Group rootGroup;
    BigDecimal epsilon;
    protected double borderThickness = 0.5;
    protected double centerX;
    protected double centerY;
    protected double innerRadius = 30;
    protected double outerRadius =60;
    protected double startAngle;
    protected double length;
    protected long seedNumber = 10;
    protected double occupancy = -1;
    protected double minLeafLower;
    protected double overallLower;
    protected double ratioToMaxLeaf;
    protected double confLowerBound;
    protected double confUpperBound;
    protected Arc innerRing;
    protected Arc outerRing;
    protected List<Arc> bands;
    protected static boolean drawTree = false;
    protected ColorStyle colorStyle = ColorStyle.occupancy;
    protected static List<Double> maxLevelOccupancies = new ArrayList<>();
    protected double entropy;
    protected static List<Double> maxLevelEntropies = new ArrayList<>();
    protected static List<Double> minLevelEntropies = new ArrayList<>();

    public enum ColorStyle {
        differenceFromEnergy,
        logOccupancy, occupancy,
        byEntropy
    }
    public void setChildren(List<KStarTreeNode> children){
        this.children = children;
    }

    public void setColorStyle(ColorStyle style) {
        colorStyle = style;
        if(children == null || children. size() < 1)
            return;
        for(KStarTreeNode child: children)
            child.setColorStyle(colorStyle);
        recolor();
    }

    private void recolor() {
        if(children == null || children. size() < 1)
            return;
        if(bands != null && bands.size() > 0){
            for(int i = 0; i < bands.size(); i++){
                bands.get(i).setFill(children.get(i).getWeightedColor());
            }
        }
    }

    public static KStarTreeNode parseTree(File file, boolean render, Map<Integer,BigDecimal> zCutoffsByLevel)
    {
        try {
            BufferedReader fileStream = new BufferedReader(new FileReader(file));
            KStarTreeNode.Builder builder = new KStarTreeNode.Builder();
            builder.setEpsilon(0.68);
            builder.setRender(render);
            AtomicLong numNodes = new AtomicLong(0);
            fileStream.lines().forEach(line -> {
                boolean addedNode = builder.addNode(line, zCutoffsByLevel);
                if (addedNode) {
                    numNodes.incrementAndGet();
                }
            });
            log("added %d nodes from the tree", numNodes.get());
            return builder.assembleTree();
        } catch(Exception e)
        {
            System.err.println("Parse tree failed: ");
            e.printStackTrace();
        }
        return null;
    }
     public static KStarTreeNode parseTree(String fileName) {
        return parseTree(new File(fileName), false, null);
     }

     public BigDecimal getLowerBound() {
        return lowerBound;
     }

     public BigDecimal getUpperBound() {
        return upperBound;
     }


    public KStarTreeNode(int level, String[] assignments, int[] confAssignments, BigDecimal lowerBound, BigDecimal upperBound,
                         double confLowerBound, double confUpperBound, double epsilon) {
        this.level = level;
        this.assignments = assignments;
        this.confAssignments = confAssignments;
        this.upperBound = upperBound;
        this.lowerBound = lowerBound;
        this.epsilon = new BigDecimal(epsilon);
        this.bandGroup = new Group();
        this.confLowerBound = confLowerBound;
        this.confUpperBound = confUpperBound;
        this.colorSeeds = new Random[assignments.length];
        if(isRoot()) {
            this.overallUpperBound = upperBound;
        }
    }

    public KStarTreeNode assign(
    	int pos, int rc, String assignment,
		BigDecimal lowerBound, BigDecimal upperBound,
		double confLowerBound, double confUpperBound
	) {

    	// make the assignments
    	String[] assignments = this.assignments.clone();
    	assignments[pos] = assignment;
    	int[] conf = this.confAssignments.clone();
    	conf[pos] = rc;

    	KStarTreeNode child = new KStarTreeNode(
    		level + 1,
			assignments,
			conf,
			lowerBound, upperBound,
			confLowerBound, confUpperBound,
			epsilon.doubleValue()
		);

    	// propagate copies of info to the new node
    	child.overallUpperBound = this.overallUpperBound;

    	// add the child node to the tree
		addChild(child);

    	return child;
	}

	public void updateBoundsFromChildren(MathContext mathContext) {

    	if (children.isEmpty()) {
    		throw new IllegalStateException("no children from which to update bounds");
		}

		BigMath lower = new BigMath(mathContext).set(0);
		BigMath upper = new BigMath(mathContext).set(0);
		confLowerBound = Double.POSITIVE_INFINITY;
		confUpperBound = Double.NEGATIVE_INFINITY;

		for (KStarTreeNode child : children) {
			lower.add(child.lowerBound);
			upper.add(child.upperBound);
			confLowerBound = Math.min(confLowerBound, child.confLowerBound);
			confUpperBound = Math.max(confUpperBound, child.confUpperBound);
		}

		lowerBound = lower.get();
		upperBound = upper.get();

		occupancy = -1;
		computeOccupancy();
	}

    public int numStatesAtLevel(int level) {
        if(this.level == level || children == null || children.size() < 1) {
            return 1;
        }
        int subTreeSum = 0;
        for(KStarTreeNode child: children) {
            subTreeSum += child.numStatesAtLevel(level);
        }
        return subTreeSum;
    }

    public BigDecimal maxWeightedErrorBound() {
        if(children == null || children.size() < 1) {
            return upperBound.subtract(lowerBound);
        }
        BigDecimal maxChildError = BigDecimal.ZERO;
        for(KStarTreeNode child: children) {
            BigDecimal childError = child.maxWeightedErrorBound();
            if(MathTools.isLessThanOrEqual(maxChildError, childError))
                maxChildError = childError;
        }
        return maxChildError;
    }

    public double maxConfErrorBound() {
        if(children == null || children.size() < 1) {
            return Math.min(0,confUpperBound)-confLowerBound;
        }
        double maxChildError = 0;
        for(KStarTreeNode child: children) {
            maxChildError = Math.max(maxChildError, child.maxConfErrorBound());
        }
        return maxChildError;
    }

    public double computeEntropy(int maxLevel) {
        if(level >= maxLevel || children == null || children.size() < 1) {
            computeOccupancy();
            if(occupancy == 0)
                return 0;
            return -1.9891/1000.0 * occupancy * Math.log(occupancy) * numStatesAtLevel(maxLevel);
        }
        double subTreeSum = 0;
        double subTreeOccupancy = 0;
        int childCount = 0;
        for(KStarTreeNode child: children) {
            subTreeSum += child.computeEntropy(maxLevel);
            subTreeOccupancy+=child.occupancy;
            if (Double.isNaN(subTreeSum)) {
                System.out.println("wtf?");
            }
        }
        return subTreeSum;
    }

    private void computeOccupancy() {
        if(occupancy < 0)
            occupancy = upperBound.divide(overallUpperBound,70, RoundingMode.HALF_UP).doubleValue();
    }

    public double computeEntropy() {
        return computeEntropy(Integer.MAX_VALUE);
    }

    public double computeEnthalpyWithEnergiesFrom(Map<String, Double> energyMap, int maxLevel) {
        if(level >= maxLevel || children == null || children.size() < 1) {
            computeOccupancy();
            String formattedAssignments = formatAssignmentsVisual(this.assignments);
            if(!energyMap.containsKey(formattedAssignments))
                return 0;
            return occupancy*energyMap.get(formattedAssignments);
        }
        double subtreeSum = 0;
        for(KStarTreeNode child: children) {
            subtreeSum+=child.computeEnthalpyWithEnergiesFrom(energyMap, maxLevel);
        }
        return subtreeSum;
    }

    public Map<String, Double> computeEnergyMap(int maxLevel) {
        Map<String, Double> outputMap = new HashMap<>();
        populateEnergyMap(outputMap, maxLevel);
        return outputMap;
    }

    private void populateEnergyMap(Map<String, Double> map, int maxLevel) {
        if(level >= maxLevel || children == null || children.size() < 1) {
            map.put(this.formatAssignmentsVisual(assignments), this.confLowerBound);
            return;
        }
        for(KStarTreeNode child: children) {
            child.populateEnergyMap(map, maxLevel);
        }
    }

    public double computeEnthalpy(int maxLevel) {
        computeOccupancy();
        if(level >= maxLevel || children == null || children.size() < 1)
            return confLowerBound*occupancy;
        double subTreeSum = 0;
        for(KStarTreeNode child: children) {
            double childEnthalpy = child.computeEnthalpy(maxLevel);
            subTreeSum += childEnthalpy;
            if(occupancy > 0 && subTreeSum < confLowerBound) {
                double subtreeOccupancy = 0;
                double subtreeEnthalpy = 0;
                for(KStarTreeNode child2: children) {
                    System.out.println(child2+":"+child2.occupancy+"->"+child2.computeEnthalpy(maxLevel));
                    subtreeEnthalpy+= child2.computeEnthalpy(maxLevel);
                    subtreeOccupancy+=child2.occupancy;
                }
                System.out.println("Subtree occupancy:"+subtreeOccupancy);
                System.out.println("Subtree enthalpy:"+subtreeEnthalpy);
                System.err.println("???");
            }
        }
        if(occupancy > 0 && subTreeSum < confLowerBound)
            System.err.println("???");
        return subTreeSum;
    }

    public double computeEnthalpy() {
        return computeEnthalpy(Integer.MAX_VALUE);
    }

    public void initStatText() {
        this.statText = new Text(toStringVisual());
        statText.setFill(Color.WHITE);
        statText.setStroke(Color.BLACK);
        statText.setFont(Font.font("Verdana", FontWeight.BOLD, 18));
        statText.setVisible(false);
    }

    public void setRender(boolean render) {
        this.drawTree = render;
    }


    public void preprocess() {
        if(level != 0){
            System.err.println("This is only run from the root for now.");
            System.exit(-1);
        }
        double minLeafLower = getMinLeafLower();
        this.ratioToMaxLeaf = 1;
        setMinLeafLower(minLeafLower);
        computeLevelMaxOccupancies();
        computeLevelMaxEntropies();
        computeLevelMinEntropies();
    }

    private void setMinLeafLower(double treeLowerBound) {
        this.overallLower = treeLowerBound;
        this.minLeafLower = getMinLeafLower();
        double colorThreshhold = 2;
        this.ratioToMaxLeaf = (colorThreshhold-Math.min(colorThreshhold,minLeafLower - overallLower))/colorThreshhold;
        if(children == null || children.size() < 1)
            return;
        for(KStarTreeNode child: children)
            child.setMinLeafLower(treeLowerBound);

    }


    private double getMinLeafLower() {
        if(children == null || children.size() < 1)
            return confLowerBound;
        double minLower = Double.POSITIVE_INFINITY;
        for(KStarTreeNode child: children) {
            double childMinLower = child.getMinLeafLower();
            if(childMinLower < minLower)
                minLower = childMinLower;
        }
        return minLower;
    }


    private Random getColorSeed() {
        if(colorSeeds[level] == null)
            colorSeeds[level] = new Random(seedNumber);
        return colorSeeds[level];
    }

    protected void setBand(double startAngle, double length) {
        this.startAngle = startAngle;
        this.length = length;
    }

    public void showRoot() {
        if(isRoot())
            this.visible = true;
    }

    private void renderCircle() {
        renderBand(centerX, centerY, innerRadius, outerRadius, 0, 360,
                computeWeights());
    }

    public void setGroup(Group g) {
        rootGroup = g;
        if(children.size() < 1)
            return;
        for(KStarTreeNode child: children)
            child.setGroup(g);

    }

    protected void renderBand(double centerX, double centerY,
                            double innerRadius, double outerRadius,
                            double arcStart, double arcLength, double[] weights) {
        if(children.size() < 1 || bandGroup.getChildren().size() > 0)
            return;
        bands = new ArrayList<>();
        ArrayList<Node> separators = new ArrayList<>();
        outerRing = new Arc(centerX, centerY, outerRadius, outerRadius, arcStart, arcLength);
        outerRing.setType(ArcType.ROUND);
        outerRing.setFill(Color.WHITE);
        innerRing = new Arc(centerX, centerY, innerRadius+borderThickness, innerRadius+borderThickness,
                arcStart, arcLength);
        innerRing.setType(ArcType.ROUND);
        innerRing.setFill(Color.WHITE);

        double startAngle = arcStart;
        final double arcLengthFinal = arcLength;// - borderThickness*(weights.length);
        for(int i = 0; i < weights.length; i++) {
            double weight = weights[i];
            double weightLength = arcLengthFinal*weight;
            Arc arc = new Arc(centerX, centerY, outerRadius-borderThickness,
                    outerRadius-borderThickness, startAngle, arcLengthFinal*weight);
            KStarTreeNode child = children.get(i);
            Color arcColor = child.getWeightedColor();
            arc.setFill(arcColor);
            arc.setType(ArcType.ROUND);
            bands.add(arc);
            double s = Math.toRadians(startAngle+weightLength);
            Line separator = new Line(centerX,//+innerRadius*Math.cos(s),
                    centerY,//-innerRadius*Math.sin(s),
                    centerX+outerRadius*Math.cos(s),
                    centerY-outerRadius*Math.sin(s));
            separator.setStroke(Color.gray(0.92));
            separator.setStrokeWidth(borderThickness);
            separators.add(separator);
            double finalStartAngle = startAngle;
            double finalLength = arcLengthFinal*weight;
            child.setBand(finalStartAngle, finalLength);
            child.innerRadius = innerRadius;
            child.outerRadius = outerRadius;
            arc.setOnMouseClicked((e)-> {
                if(e.isControlDown())
                    child.toggleBand();
            });
            arc.setOnMouseEntered((e)-> {
                arc.setFill(child.getWeightedColor().brighter()
                        .desaturate()
                        .desaturate());
                child.showConfInfo(e.getX()+10, e.getY());
            });
            arc.setOnMouseMoved((e) -> {
                child.showConfInfo(e.getX()+10, e.getY());
            });
            arc.setOnMouseExited((e)-> {
                arc.setFill(child.getWeightedColor());
                child.hideConfInfo();
            });
            startAngle += weight*arcLengthFinal;
        }
        bandGroup.getChildren().addAll(outerRing, innerRing);
        bandGroup.getChildren().addAll(bands);
        bandGroup.getChildren().addAll(separators);
        bandGroup.setVisible(true);
        rootGroup.getChildren().add(bandGroup);
        innerRing.toBack();
        for(Node sep:separators)
            sep.toBack();
        for(Node band: bands)
            band.toBack();
        outerRing.toBack();
        bandGroup.toBack();

    }

    public void computeLevelMaxOccupancies() {
        if(children == null || children.size() < 1)
            return;
        while(maxLevelOccupancies.size() <= level+1)
            maxLevelOccupancies.add(Double.NEGATIVE_INFINITY);
        for(KStarTreeNode child: children) {
            maxLevelOccupancies.set(level+1, Math.max(child.occupancy, maxLevelOccupancies.get(level+1)));
            child.computeLevelMaxOccupancies();
        }

    }

    double energyThreshold = 0.5;
    double ratioThreshold = 0.8;
    double redblueEnergyThreshold = 2;
    double occupancyThreshold = 0.5;
    double logOccupancyThreshold = 0.80;
    double entropyThreshold = 0.5;

    protected Color getWeightedColor() {
        switch(colorStyle) {
            case differenceFromEnergy:
                return getEnergyWeightedColor();
            case occupancy:
                return getOccupancyWeightedColor();
            case logOccupancy:
                return getLogOccupancyWeightedColor();
            case byEntropy:
                return getEntropyWeightedColor();
            default:
                return getOccupancyWeightedColor();

        }
    }

    public int numConfsWithin(double diffFromGMEC) {
        if(children == null || children.size()<1
                && confLowerBound - overallLower < diffFromGMEC) {
            return 1;
        }
        int sum = 0;
        for(KStarTreeNode child: children)
            sum+=child.numConfsWithin(diffFromGMEC);
        return sum;

    }

    public double[] computeEnergyErrorWithinEnergyRange(double diffFromGMEC) {
        double[] range = new double[2];
        computeEnergyErrorWithinEnergyRange(diffFromGMEC, range);
        return range;
    }

    private void computeEnergyErrorWithinEnergyRange(double diffFromGMEC, double[] range) {
        if(children == null || children.size()<1
            && confLowerBound - overallLower < diffFromGMEC) {
            double diff = range[1] - range[0];
            if (confUpperBound - confLowerBound > diff || range[0] == 0) {
                range[0] = confLowerBound;
                range[1] = confUpperBound;
            }
        }
        for(KStarTreeNode child: children)
            child.computeEnergyErrorWithinEnergyRange(diffFromGMEC, range);

    }

    protected Color getLogOccupancyWeightedColor() {
        double logOccupancy = Math.log(occupancy);
        double occupancy = Math.max(0.00000000000000000000000001,1-logOccupancy/Math.log(0.001));
        if(occupancy < logOccupancyThreshold)
            return redBlueGradient(occupancy/logOccupancyThreshold);
        double weight = (occupancy - logOccupancyThreshold)/(1-logOccupancyThreshold);
        return blueGreenGradient(weight);
    }


    protected Color getOccupancyWeightedColor() {
        double levelMaxOccupancy = maxLevelOccupancies.get(level);
        double occupancy = this.occupancy/levelMaxOccupancy;
        if(occupancy < occupancyThreshold)
            return redBlueGradient(occupancy/occupancyThreshold);
        double weight = (occupancy - occupancyThreshold)/(1-occupancyThreshold);
        return blueGreenGradient(weight);
    }

    protected Color getEnergyWeightedColor() {
        if(minLeafLower - overallLower > energyThreshold)
            return redBlueGradient((redblueEnergyThreshold-Math.min(minLeafLower-overallLower, redblueEnergyThreshold))/redblueEnergyThreshold);
        double energyWeight = (minLeafLower -overallLower)/0.5;
        return blueGreenGradient(1-energyWeight);
    }


    private Color redGreenGradient(double ratioToMaxLeaf) {
        double newRatio = ratioToMaxLeaf;
        return new Color(1-newRatio, newRatio, 0, 1);
    }

    protected Color blueGreenGradient(double ratioToMaxLeaf) {
        double newRatio = ratioToMaxLeaf;
        return new Color(0.15*newRatio, newRatio, 0.5*(1-newRatio)+0.2, 1);
    }

    protected Color redBlueGradient(double ratioToMaxLeaf) {
        double newRatio = ratioToMaxLeaf;
        return new Color(1*(1-newRatio), 0, newRatio,1);
    }

    protected void hideConfInfo() {
        statText.setVisible(false);
    }

    protected void showConfInfo(double mouseX, double mouseY) {
        statText.setTranslateX(mouseX);
        statText.setTranslateY(mouseY);
        statText.setText(this.toStringVisual());
        statText.setVisible(true);
    }

    protected void toggleBand() {
        expanded = !expanded;
        if(hasChildren() && !childrenRendered){
            expanded = true;
            renderBand(centerX, centerY, outerRadius-borderThickness, outerRadius+20,
                    startAngle, length, computeWeights());
            childrenRendered = true;
            expanded = true;
            setVisible(true);
        }
        setVisible(expanded);
    }


    protected double[] computeWeights() {
        return normalizeWeights();
    }

    private double[] normalizeWeights() {
        if(children == null || children.size() < 1)
            return new double[]{1.0};
        double[] weights = new double[children.size()];
        IntStream.range(0,children.size()).forEach(i->
                weights[i] = children.get(i).upperBound.divide(upperBound,5,RoundingMode.HALF_UP).doubleValue()
        );
        return weights;
    }

    private double[] threshold(double[] weights, double v) {
        int numAbovethreshold = 0;
        for(int i =0; i < weights.length; i++) {
            if(weights[i] > v)
                numAbovethreshold++;
        }

        double[] out = new double[numAbovethreshold];
        double leftover = 1;
        for(int i = 0; i<out.length;i++) {
            out[i] = weights[i];
            leftover-=weights[i];
        }

        return out;
    }


    private double[] toPrimitive(ObservableList<Double> source)
    {
        double[] output = new double[source.size()];
        IntStream.range(0, source.size())
                .forEach(i -> output[i] = (double) source.get(i));
        return output;
    }

    private Double[] toObject(double[] source)
    {
        Double[] output = new Double[source.length];
        IntStream.range(0, source.length)
                .forEach(i -> output[i] = (Double) source[i]);
        return output;
    }

    public void recenter(double x, double y) {
        centerX = x;
        centerY = y;
        bandGroup.relocate(x,y);

        redraw();
    }

    private void redraw() {
    }

    public void autoExpand(double v) {
        autoExpand(v, Integer.MAX_VALUE);
    }

    public void autoExpand(double v, int maxLevel) {
        if(level >= maxLevel) {
            setChildrenVisible(false);
            return;
        }
        if(occupancy > v) {
            toggleBand();
        }
        if(children == null || children.size() < 1)
            return;
        for(KStarTreeNode child: children) {
            child.autoExpand(v, maxLevel);
        }
    }



    public Map<KStarTreeNode,List<KStarTreeNode>> getTopSamples(int numSamples, int levelThreshold) {
        Map<KStarTreeNode,List<KStarTreeNode>> samples = new HashMap<>();
        getTopSamples(numSamples, levelThreshold, samples);
        return samples;
    }

    public List<KStarTreeNode> getTopSamplesInSubtree(int numSubtreeSamples) {
        PriorityQueue<KStarTreeNode> queue = new PriorityQueue<>();
        List<KStarTreeNode> confs = new ArrayList<>();
        queue.add(this);
        while(!queue.isEmpty() && confs.size() < numSubtreeSamples) {
            KStarTreeNode node = queue.poll();
            if(node.fullyAssigned())
                confs.add(node);
            else {
                queue.addAll(node.children);
            }
        }
        return confs;
    }

    public void getTopSamples(int numSamples, int levelThreshold, Map<KStarTreeNode, List<KStarTreeNode>> lists) {
        System.out.println("At "+this+", looking for "+numSamples+" samples");
        if(numSamples < 1) {
            System.out.println("No samples needed here.");
            return;
        }
        if(children == null || children.size() < 1) {
            System.out.println("No children. Done here.");
            if(fullyAssigned()) {
                List<KStarTreeNode> thisList = new ArrayList<>();
                thisList.add(this);
                lists.put(this, thisList);
            }
            return;
        }
        if(this.level >= levelThreshold) {
            List<KStarTreeNode> subtreeList = new ArrayList<>();
            lists.put(this,getTopSamplesInSubtree(numSamples));
            return;
        }


        int numNonZeroChildren = 0;
        for(KStarTreeNode child: children) {
            double childOccupancy = child.getOccupancy();
            if (childOccupancy > 0.01)
                numNonZeroChildren++;
        }
        if(numNonZeroChildren < 1)
            return;
        int samplesPerChild = numSamples/numNonZeroChildren;
        int leftovers = numSamples%numNonZeroChildren;
        for(int i = 0; i < numNonZeroChildren; i++) {
            if(i < leftovers)
                children.get(i).getTopSamples(samplesPerChild+1, levelThreshold, lists);
            else {
                children.get(i).getTopSamples(samplesPerChild, levelThreshold, lists);
                if(samplesPerChild < 1)
                    return;
            }
        }

    }

    private double getOccupancy() {
        if(occupancy < 0)
            computeOccupancy();
        return occupancy;
    }

    private boolean fullyAssigned() {
        boolean fullyAssigned = true;
        for(String assignment: assignments) {
            if(assignment.contains("*"))
                fullyAssigned = false;
        }
        return fullyAssigned;
    }

    public String[] getAssignments() {
        return assignments;
    }

    public int[] getConfAssignments() {
        return confAssignments;
    }

    public double getConfLowerBound() {
        return confLowerBound;
    }
    public double getConfUpperBound() {
        return confUpperBound;
    }

    public String getEnsemblePDBName() {
        String out = "";
        for(String i: assignments) {
            i = i.replace(":", "-");
            if (!i.contains("*"))
                out += i + "-";
        }
        return out;
    }

    public void pieChart(int... levels) {
        Set<Integer> levelSet = new HashSet<>();
        for(int level: levels)
            levelSet.add(level);
        setVisibleLevels(levelSet);
    }

    private void setVisibleLevels(Set<Integer> levelSet) {
        if(!levelSet.contains(level+1)) {
            if(innerRing != null)
                innerRing.setVisible(false);
            if(outerRing != null)
                outerRing.setVisible(false);
            if(bands != null) {
               for(Node band: bands)
                   band.setVisible(false);
            }
        }
        if(levelSet.contains(level +1)) {
            this.visible = true;
            if(bands != null) {
                for(Node band: bands)
                    band.setVisible(true);
            }
            if(innerRing!=null)
                innerRing.setVisible(false);
            if(outerRing != null)
                outerRing.setVisible(true);
        }

        if(children == null || children.size() <1 )
            return;
        for(KStarTreeNode child: children)
            child.setVisibleLevels(levelSet);
    }

    public void pieChart(int targetLevel) {
        if(level +1 < targetLevel) {
            if(innerRing != null)
                innerRing.setVisible(false);
            if(outerRing != null)
                outerRing.setVisible(false);
            if(bands != null) {
               for(Node band: bands)
                   band.setVisible(false);
            }
        }
        if(level +1 == targetLevel) {
            this.visible = true;
            if(bands != null) {
                for(Node band: bands)
                    band.setVisible(true);
            }
            if(innerRing!=null)
                innerRing.setVisible(true);
            if(outerRing != null)
                outerRing.setVisible(true);
        }
        if(level +1 > targetLevel) {
            if(bands != null) {
                for(Node band: bands)
                    band.setVisible(false);
            }
            if(innerRing!=null)
                innerRing.setVisible(false);
            if(outerRing != null)
                outerRing.setVisible(false);
        }

        if(children == null || children.size() <1 )
            return;
        for(KStarTreeNode child: children)
            child.pieChart(targetLevel);
    }

    public void setTextRoot(Group textGroup) {
        initStatText();
        textGroup.getChildren().add(statText);
        textGroup.setVisible(true);
    }

    public void toggleCenter() {
        if(innerRing!=null)
            innerRing.setVisible(!innerRing.isVisible());
    }

    public Set<KStarTreeNode> getLevelNodes(int targetLevel) {
        HashSet<KStarTreeNode> nodes = new HashSet<>();
        getLevelNodes(targetLevel, nodes);
        return nodes;
    }

    private void getLevelNodes(int targetLevel, Set<KStarTreeNode> nodes) {
        if(this.level == targetLevel)
            nodes.add(this);
        if(level > targetLevel || children == null || children.size() < 1)
            return;
        for(KStarTreeNode child: children)
            child.getLevelNodes(targetLevel, nodes);
    }


    public static class Builder
    {
        protected KStarTreeNode root;
        protected Stack<KStarTreeNode> buildStack = new Stack<>();
        protected int lastLevel = -1;
        protected KStarTreeNode lastNode;
        protected double epsilon = 0;
        protected boolean render = false;

        public Builder setRender(boolean render) {
            this.render = render;
            return this;
        }
        public boolean addNode(String line, Map<Integer,BigDecimal> zCutoffsByLevel)
        {
            // regex all the things!!
            Matcher m = p.matcher(line);
            if(m.matches() && debug) {
                System.out.println("Groups:");
                for(int i = 0; i < m.groupCount(); i++)
                    System.out.println("group "+i+":"+m.group(i));
            }
            if(!m.matches())
                System.out.println(line);
            int level = m.group(1).length()/2;
            String[] bounds = m.group(6).split(",");
            int[] confAssignments = Arrays.stream(m.group(3).replaceAll(" ","").split(",")).mapToInt(Integer::parseInt).toArray();
            String[] assignments = m.group(4).split(",");
            BigDecimal lowerBound = new BigDecimal(bounds[0]);
            BigDecimal upperBound = new BigDecimal(bounds[1]);
            String[] confBounds = m.group(5).split(",");
            double confLowerBound = Double.valueOf(confBounds[0]);
            double confUpperBound = Double.valueOf(confBounds[1]);

            if (zCutoffsByLevel != null) {

                // if this node is too "small", just drop it entirely
                if (upperBound.compareTo(zCutoffsByLevel.get(level)) < 0) {
                    return false;
                }

                // if this node had its parent dropped, drop it too
                if (level > lastLevel + 1) {
                    return false;
                }
            }

            if(level > lastLevel) {
                buildStack.push(lastNode);
            }
            while(level < lastLevel) {
                lastLevel--;
                buildStack.pop();
            }
            KStarTreeNode newNode = new KStarTreeNode(level, assignments, confAssignments, lowerBound, upperBound,
                    confLowerBound, confUpperBound, epsilon);
            KStarTreeNode curParent = buildStack.peek();
            if(newNode.isRoot()) {
                root = newNode;
                if(render)
                    root.initStatText();
            }
            else
                curParent.addChild(newNode);
            lastLevel = level;
            lastNode = newNode;

            return true;
        }

        public Builder setEpsilon(double epsilon)
        {
            this.epsilon = epsilon;
            return this;
        }

        public KStarTreeNode assembleTree()
        {
            root.prepTree();
            return root;
        }


    }

    protected void prepTree()
    {
        if(children == null)
            return;
        if(children != null && children.size() > 0)
        {
            Collections.sort(children);
            for(KStarTreeNode child: children)
                child.prepTree();
        }
    }

    public boolean isChildOf(KStarTreeNode otherNode) {
        if(parent == null)
            return false;
        if(parent == otherNode)
            return true;
        return parent.isChildOf(otherNode);
    }

    public boolean isParentOf(KStarTreeNode otherNode) {
        return otherNode.isChildOf(this);
    }

    protected void addChild(KStarTreeNode newNode) {

    	if (newNode.parent != null) {
    		throw new IllegalArgumentException("node already has a parent");
		}

        children.add(newNode);
        newNode.statText = this.statText;
        newNode.parent = this;
        newNode.overallUpperBound = this.overallUpperBound;
        newNode.computeOccupancy();
    }

    public KStarTreeNode parent() {
    	return parent;
	}

    public void removeFromParent() {
    	if (parent != null) {
    		boolean wasRemoved = parent.children.remove(this);
    		if (!wasRemoved) {
    			throw new IllegalStateException("node was not a child of its parent");
			}
    		parent = null;
		}
	}

	public Integer getAssignmentIndex() {

    	// root node? no assignment
		if (parent == null) {
			return null;
		}

		// we should have exactly one more assignment than our parent
		// which one is it?
		for (int i=0; i<assignments.length; i++) {
			int myrc = confAssignments[i];
			int parentrc = parent.confAssignments[i];
			if (myrc != Conf.Unassigned && parentrc == Conf.Unassigned) {
				return i;
			}
		}

		// something is wrong, find a programmer
		throw new Error("tree structure doesn't match assignment order. this is a bug.");
	}

    public String toString()
    {
        return formatAssignments(assignments)+":["+confLowerBound+","+confUpperBound+"]->["+formatBigDecimal(lowerBound)+", "+formatBigDecimal(upperBound)+"]";
    }

    public String toStringVisual()
    {
        return formatAssignmentsVisual(assignments)+":"+formatBigDecimal(lowerBound)+", "+formatBigDecimal(upperBound);
    }


    private String formatBigDecimal(BigDecimal bd) {
        return format(bd);
    }

    private String formatAssignmentsVisual(String[] assignments) {
        String out = "[";
        for(String i: assignments)
            if(!i.contains("*"))
                out+=i+",\n";
        return out+"]";
    }

    private String formatAssignments(String[] assignments) {
        String out = "[";
        for(String i: assignments)
            if(!i.contains("*"))
                out+=i+",";
            else out+="*,";
        return out+"]";
    }

    private BigInteger numSignificantConfs() {
        if(children == null || children.size() < 1) {

            if (upperBound.compareTo(BigDecimal.ZERO) > 0 &&
                    upperBound.divide(overallUpperBound,10, RoundingMode.HALF_UP)
                            .compareTo(epsilon)>0) {
                return BigInteger.ONE;
            } else {
                return BigInteger.ZERO;
            }
        }
        BigInteger sum = BigInteger.ZERO;
        for (KStarTreeNode child : children)
            sum = sum.add(child.numSignificantConfs());
        return sum;
    }

    protected boolean isRoot() {
        return level == 0;
    }

    private String toRenderString()
    {
        String output = formatAssignments(this.assignments);
        if(hasChildren())
            output+=numSignificantConfs().toString()+":";
        output += formatBigDecimal(lowerBound)+","+formatBigDecimal(upperBound);
        return output;

    }

    public void render() {
        if(isRoot())
            renderCircle();
    }
    private static String format(BigDecimal x)
    {
        NumberFormat formatter = new DecimalFormat("0.0E0");
        formatter.setRoundingMode(RoundingMode.HALF_UP);
        //formatter.setMinimumFractionDigits((x.scale() > 0) ? x.precision() : x.scale());
        formatter.setMaximumFractionDigits(3);
        return formatter.format(x);
    }

    private void setVisible(boolean visible) {
        if(this.visible != visible) {
            this.visible = visible;
            if(innerRing != null)
                innerRing.setVisible(visible);
            if(bands != null)
                for(Node band: bands)
                    band.setVisible(visible);
        }
        if(bandGroup!=null)
            bandGroup.setVisible(visible);
        if(hasChildren())
            for(KStarTreeNode child : children)
                child.setVisible(visible);
    }

    @Override
    public int compareTo(KStarTreeNode other)
    {
        return -upperBound.compareTo(other.upperBound);
    }

    private boolean hasChildren() {
        return children != null && children.size() > 0;
    }

    private void setChildrenVisible(boolean visible) {
        expanded = visible;
        if(hasChildren())
            for(KStarTreeNode child : children)
                child.setVisible(visible);
    }

    public void render(Group g)
    {
        render();
    }

    public void printTree(String prefix, FileWriter writer)
    {
        String confString = toString();
        String out = prefix+confString;
        if(writer != null) {
            try {
                writer.write(out+"\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        else
            System.out.println(out);
        if(children != null && !children.isEmpty()) {
            Collections.sort(children, (a,b)-> -a.upperBound
                    .compareTo(b.upperBound));
            for (KStarTreeNode child : children)
                child.printTree(prefix + "~+", writer);
        }
    }

    public void printTree() {
        printTree("", null);
    }

    /** printTree uses a different format than MARK* apparently */
    public void printTreeLikeMARKStar(Writer out)
	throws IOException {
    	printTreeLikeMARKStar(out, "");
	}

	public void printTreeLikeMARKStar(Writer out, String prefix)
	throws IOException {

    	// render the prefix
		out.append(prefix);

    	// render the conf
		out.append("(");
		for (int rc : confAssignments) {
			out.append(Integer.toString(rc));
			out.append(", ");
		}
		out.append(")");

		out.append("->");

		// render the assignments
		out.append("(");
		for (int i=0; i<assignments.length; i++) {
			if (i > 0) {
				out.append(",");
			}
			out.append(assignments[i]);
		}
		out.append(")");

		out.append(":");

		// render the conf bounds
		out.append("[");
		out.append(Double.toString(confLowerBound));
		out.append(",");
		out.append(Double.toString(confUpperBound));
		out.append("]");

		out.append("->");

		// render the zSum bounds
		out.append("[");
		out.append(lowerBound.toString()); // this will print an exact representation of the BigDecimal
		out.append(",");
		out.append(upperBound.toString());
		out.append("]");

		out.append("\n");

		// recurse
		for (KStarTreeNode child : children) {
			child.printTreeLikeMARKStar(out, prefix + "~+");
		}
	}

    protected Color getEntropyWeightedColor() {
        double levelMaxEntropy = maxLevelEntropies.get(level);
        double levelMinEntropy = minLevelEntropies.get(level);
        double ent = (this.entropy - levelMinEntropy) / (levelMaxEntropy - levelMinEntropy);
        if(ent < entropyThreshold)
            return redBlueGradient(ent/entropyThreshold);
        double weight = (ent - entropyThreshold)/(1-entropyThreshold);
        return blueGreenGradient(weight);
    }

    public void computeLevelMaxEntropies() {
        if(children == null || children.size() < 1)
            return;
        while(maxLevelEntropies.size() <= level+1)
            maxLevelEntropies.add(Double.NEGATIVE_INFINITY);
        for(KStarTreeNode child: children) {
            SeqTreeNode seqNode = (SeqTreeNode) child;
            maxLevelEntropies.set(level+1, Math.max((seqNode).getEntropy(), maxLevelEntropies.get(level+1)));
            seqNode.computeLevelMaxEntropies();
        }

    }

    public void computeLevelMinEntropies() {
        if(children == null || children.size() < 1)
            return;
        while(minLevelEntropies.size() <= level+1)
            minLevelEntropies.add(Double.POSITIVE_INFINITY);
        for(KStarTreeNode child: children) {
            SeqTreeNode seqNode = (SeqTreeNode) child;
            if(seqNode.getEntropy() > 0)
                minLevelEntropies.set(level+1, Math.min(seqNode.getEntropy(), minLevelEntropies.get(level+1)));
            seqNode.computeLevelMinEntropies();
        }

    }

    public void setEntropy(double val){
        this.entropy = val;
    }

    public double getEntropy(){
        return this.entropy;
    }
}
