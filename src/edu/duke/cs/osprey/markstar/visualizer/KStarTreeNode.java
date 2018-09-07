package edu.duke.cs.osprey.markstar.visualizer;

import javafx.collections.ObservableList;
import javafx.geometry.Point2D;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.paint.Color;
import javafx.scene.shape.Arc;
import javafx.scene.shape.ArcType;
import javafx.scene.shape.Line;
import javafx.scene.shape.Polygon;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.scene.transform.Transform;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.IntStream;

public class KStarTreeNode implements Comparable<KStarTreeNode>{
    public static final Pattern p = Pattern.compile("((~\\+)*)\\((.*)\\)->\\((.*)\\)\\: ?\\[(.*)\\]->\\[(.*)\\](.*)");
    private static final boolean debug = false;
    private static Random[] colorSeeds;
    private BigDecimal overallUpperBound;
    private int level = -1;
    private String[] assignments;
    private int[] confAssignments;
    private BigDecimal upperBound;
    private BigDecimal lowerBound;
    private List<KStarTreeNode> children = new ArrayList<>();
    private KStarTreeNode parent;
    private boolean visible =  false;
    private boolean expanded = false;
    private boolean childrenRendered = false;
    private Text statText;
    private Group bandGroup;
    private Group rootGroup;
    private BigDecimal epsilon;
    private double borderThickness = 0.5;
    private double centerX;
    private double centerY;
    private double innerRadius = 30;
    private double outerRadius =60;
    private double startAngle;
    private double length;
    private long seedNumber = 10;
    private double occupancy = 0;
    private double minLeafLower;
    private double overallLower;
    private double ratioToMaxLeaf;
    private double confLowerBound;
    private double confUpperBound;
    private static boolean drawTree = false;

    public static KStarTreeNode parseTree(File file, boolean render)
    {
        try {
            BufferedReader fileStream = new BufferedReader(new FileReader(file));
            KStarTreeNode.Builder builder = new KStarTreeNode.Builder();
            builder.setEpsilon(0.68);
            builder.setRender(render);
            fileStream.lines().forEach(line -> builder.addNode(line));
            return builder.assembleTree();
        } catch(Exception e)
        {
            System.err.println("Parse tree failed: ");
            e.printStackTrace();
        }
        return null;
    }
     public static KStarTreeNode parseTree(String fileName) {
        return parseTree(new File(fileName), false);
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

    public void initStatText() {
        this.statText = new Text(toStringVisual());
        statText.setFill(Color.WHITE);
        statText.setStroke(Color.BLACK);
        statText.setFont(Font.font("Verdana", FontWeight.BOLD, 12));
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

    private void setBand(double startAngle, double length) {
        this.startAngle = startAngle;
        this.length = length;
    }

    public void showRoot() {
        if(isRoot())
            this.visible = true;
    }

    private void renderCircle(Group g) {
        renderBand(g, centerX, centerY, innerRadius, outerRadius, 0, 360,
                computeWeights());
    }

    public void setGroup(Group g) {
        rootGroup = g;
        if(children.size() < 1)
            return;
        for(KStarTreeNode child: children)
            child.setGroup(g);

    }

    private void renderBand(Group g, double centerX, double centerY,
                            double innerRadius, double outerRadius,
                            double arcStart, double arcLength, double[] weights) {
        if(children.size() < 1 || bandGroup.getChildren().size() > 0)
            return;
        ArrayList<Node> bands = new ArrayList<>();
        ArrayList<Node> separators = new ArrayList<>();
        Arc backRing = new Arc(centerX, centerY, outerRadius, outerRadius, arcStart, arcLength);
        backRing.setType(ArcType.ROUND);
        backRing.setFill(Color.WHITE);
        Arc innerRing = new Arc(centerX, centerY, innerRadius+borderThickness, innerRadius+borderThickness,
                arcStart, arcLength);
        innerRing.setType(ArcType.ROUND);
        innerRing.setFill(Color.WHITE);

        Arc whiteCenter = new Arc(centerX, centerY, innerRadius, innerRadius, arcStart, arcLength);
        whiteCenter.setType(ArcType.ROUND);
        whiteCenter.setFill(Color.WHITE);
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
            Line separator = new Line(centerX+innerRadius*Math.cos(s),
                    centerY-innerRadius*Math.sin(s),
                    centerX+outerRadius*Math.cos(s),
                    centerY-outerRadius*Math.sin(s));
            separator.setStroke(Color.WHITE);
            separator.setStrokeWidth(borderThickness);
            separators.add(separator);
            double finalStartAngle = startAngle;
            double finalLength = arcLengthFinal*weight;
            child.setBand(finalStartAngle, finalLength);
            child.innerRadius = innerRadius;
            child.outerRadius = outerRadius;
            arc.setOnMouseClicked((e)-> {
                child.toggleBand();
                System.out.println("Clicked "+arc+", visible? :"+child.visible);
            });
            arc.setOnMouseEntered((e)-> {
                arc.setFill(arcColor.brighter()
                        .desaturate()
                        .desaturate());
                child.showConfInfo(e.getX()+10, e.getY());
            });
            arc.setOnMouseMoved((e) -> {
                child.showConfInfo(e.getX()+10, e.getY());
            });
            arc.setOnMouseExited((e)-> {
                arc.setFill(arcColor);
                child.hideConfInfo();
            });
            startAngle += weight*arcLengthFinal;
        }
        bandGroup.getChildren().addAll(backRing, innerRing, whiteCenter);
        bandGroup.getChildren().addAll(bands);
        bandGroup.getChildren().addAll(separators);
        bandGroup.setVisible(true);
        g.getChildren().add(bandGroup);
        whiteCenter.toBack();
        innerRing.toBack();
        for(Node sep:separators)
            sep.toBack();
        for(Node band: bands)
            band.toBack();
        backRing.toBack();
        bandGroup.toBack();

    }

    double energyThreshold = 0.5;
    double ratioThreshold = 0.8;
    double redblueEnergyThreshold = 2;

    private Color getWeightedColor() {
        if(minLeafLower - overallLower > energyThreshold)
            return redBlueGradient((redblueEnergyThreshold-Math.min(minLeafLower-overallLower, redblueEnergyThreshold))/redblueEnergyThreshold);
        double energyWeight = (minLeafLower -overallLower)/0.5;
        return blueGreenGradient(1-energyWeight);
    }

    private Color redGreenGradient(double ratioToMaxLeaf) {
        double newRatio = ratioToMaxLeaf;
        return new Color(1-newRatio, newRatio, 0, 1);
    }

    private Color blueGreenGradient(double ratioToMaxLeaf) {
        double newRatio = ratioToMaxLeaf;
        return new Color(0.15*newRatio, newRatio, 0.5*(1-newRatio)+0.2, 1);
    }

    private Color redBlueGradient(double ratioToMaxLeaf) {
        double newRatio = ratioToMaxLeaf;
        return new Color(1*(1-newRatio), 0, newRatio,1);
    }

    private Color getRandomColor() {
        Random colorSeed = getColorSeed();
        double baseBrightness = 0.2;
        return new Color(colorSeed.nextDouble()*(1-baseBrightness)+baseBrightness,
                colorSeed.nextDouble()*(1-baseBrightness)+baseBrightness,
                colorSeed.nextDouble()*(1-baseBrightness)+baseBrightness, 1).saturate().saturate();
    }

    private void hideConfInfo() {
        statText.setVisible(false);
    }

    private void showConfInfo(double mouseX, double mouseY) {
        statText.setX(mouseX);
        statText.setY(mouseY);
        statText.setText(this.toStringVisual());
        statText.setVisible(true);
    }

    private void toggleBand() {
        expanded = !expanded;
        if(hasChildren() && !childrenRendered){
            expanded = true;
            renderBand(rootGroup, centerX, centerY, outerRadius-borderThickness, outerRadius+20,
                    startAngle, length, computeWeights());
            childrenRendered = true;
            expanded = true;
            setVisible(true);
        }
        setVisible(expanded);
    }


    private double[] computeWeights() {
        return normalizeWeights();
    }

    private double[] normalizeWeights() {
        if(children == null || children.size() < 1)
            return new double[]{1.0};
        double[] weights = new double[children.size()];
        IntStream.range(0,children.size()).forEach(i->
                weights[i] = children.get(i).upperBound.divide(upperBound,5,RoundingMode.HALF_UP).doubleValue()
        );
        return threshold(weights,0.01);
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

    private void transformPolygon(Polygon p, Transform t)
    {
        double[] points = toPrimitive(p.getPoints());
        double[] transformed = new double[points.length];
        t.transform2DPoints(points,0,transformed,0, points.length/2);
        p.getPoints().setAll(toObject(transformed));
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
        if(!rootGroup.getChildren().contains(statText))
            rootGroup.getChildren().add(statText);
        if(occupancy > v) {
            toggleBand();
        }
        if(children == null || children.size() < 1)
            return;
        for(KStarTreeNode child: children) {
            child.autoExpand(v);
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
            double childOccupancy = child.upperBound.divide(overallUpperBound,5,RoundingMode.HALF_UP).doubleValue();
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

    public double getConfLowerbound() {
        return confLowerBound;
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

    public static class Builder
    {
        private KStarTreeNode root;
        private Stack<KStarTreeNode> buildStack = new Stack<>();
        private int lastLevel = -1;
        private KStarTreeNode lastNode;
        private double epsilon = 0;
        private boolean render = false;

        public Builder setRender(boolean render) {
            this.render = render;
            return this;
        }
        public void addNode(String line)
        {
            Matcher m = p.matcher(line);
            if(m.matches() && debug) {
                System.out.println("Groups:");
                for(int i = 0; i < m.groupCount(); i++)
                    System.out.println("group "+i+":"+m.group(i));
            }
            int level = m.group(1).length()/2;
            String[] bounds = m.group(6).split(",");
            int[] confAssignments = Arrays.stream(m.group(3).replaceAll(" ","").split(",")).mapToInt(Integer::parseInt).toArray();
            String[] assignments = m.group(4).split(",");
            BigDecimal lowerBound = new BigDecimal(bounds[0]);
            BigDecimal upperBound = new BigDecimal(bounds[1]);
            String[] confBounds = m.group(5).split(",");
            double confLowerBound = Double.valueOf(confBounds[0]);
            double confUpperBound = Double.valueOf(confBounds[1]);
            if(level > lastLevel) {
                buildStack.push(lastNode);
            }
            if(level < lastLevel) {
                while(level < lastLevel) {
                    lastLevel--;
                    buildStack.pop();
                }
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

    private void prepTree()
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

    private void addChild(KStarTreeNode newNode) {
        children.add(newNode);
        newNode.statText = this.statText;
        newNode.parent = this;
        newNode.overallUpperBound = this.overallUpperBound;
        newNode.occupancy = upperBound.divide(overallUpperBound,5,RoundingMode.HALF_UP).doubleValue();
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

    private boolean isRoot() {
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

    public void render(Group g, Point2D start) {
        if(isRoot())
            renderCircle(g);
    }
    private static String format(BigDecimal x)
    {
        NumberFormat formatter = new DecimalFormat("0.0E0");
        formatter.setRoundingMode(RoundingMode.HALF_UP);
        formatter.setMinimumFractionDigits((x.scale() > 0) ? x.precision() : x.scale());
        return formatter.format(x);
    }

    private void setVisible(boolean visible) {
        if(this.visible != visible) {
            this.visible = visible;
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
        render(g, new Point2D(0,0));
    }
}
