package edu.duke.cs.osprey.markstar.visualizer;

import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.geometry.Point2D;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import org.jetbrains.annotations.NotNull;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Arrays;
import java.util.Optional;

public class Visualizer extends Application {

    private static final double SCALE_DELTA = 1.1;
    Stage primaryStage;
    BorderPane triroot;
    Pane ringNode;
    KStarTreeNode root;
    Circle dot = new Circle(0,0,10);
    Group rootGroup;
    private double mouseDownX;
    private double mouseDownY;
    private double ringX;
    private double ringY;

    @Override
    public void start(Stage primaryStage) throws Exception{
        this.primaryStage = primaryStage;
        primaryStage.setTitle("MARK* Tree Analyzer (v0.2)");
        triroot = new BorderPane();
        triroot.setBackground(new Background(new BackgroundFill(Color.color(1,1,1), CornerRadii.EMPTY, Insets.EMPTY)));
        Scene test = new Scene(triroot, 300, 275);
        primaryStage.setScene(test);
        MenuBar menuBar = getMenuBar(primaryStage);
        triroot.setTop(menuBar);
        primaryStage.setScene(test);
        primaryStage.show();

    }

    @NotNull
    private MenuBar getMenuBar(Stage primaryStage) {
        final Menu file = new Menu("File");
        final Menu options = new Menu("Options");
        final Menu help = new Menu("Help");
        MenuItem loadTree = new MenuItem("Load tree!");
        loadTree.setOnAction(e -> {
            FileChooser fc = new FileChooser();
            File selectedFile = fc.showOpenDialog(primaryStage);
            loadTreeFromFile(selectedFile);
        });
        MenuItem helpDevShortCut = new MenuItem("DevShortCut");
        helpDevShortCut.setOnAction(e -> {
            devShortCut();
        });
        help.getItems().add(helpDevShortCut);
        MenuItem helpDevShortCut2 = new MenuItem("DevShortCut2");
        helpDevShortCut2.setOnAction(e -> {
            devShortCut2();
        });
        help.getItems().add(helpDevShortCut2);
        MenuItem helpDevShortCut3 = new MenuItem("DevShortCut3");
        helpDevShortCut3.setOnAction(e -> {
            devShortCut3();
        });
        MenuItem setvisibleLevels = new MenuItem("Set visible levels");
        setvisibleLevels.setOnAction(e-> {
            TextInputDialog dialog = new TextInputDialog("1");
            dialog.setHeaderText("Input levels to show");

            // Traditional way to get the response value.
            Optional<String> result = dialog.showAndWait();
            result.ifPresent(levelsString -> {
                int[] levels = Arrays.stream(levelsString.split(",")).mapToInt(Integer::parseInt).toArray();
                root.pieChart(levels);
            });
        });
        MenuItem colorByOccupancy = new MenuItem("Color by statistical weight");
        colorByOccupancy.setOnAction((e)->{root.setColorStyle(KStarTreeNode.ColorStyle.occupancy);});
        MenuItem colorByEnergy = new MenuItem("Color by energy");
        colorByEnergy.setOnAction((e)->{root.setColorStyle(KStarTreeNode.ColorStyle.differenceFromEnergy);});
        MenuItem colorByLogOccupancy = new MenuItem("Color by log statistical weight");
        colorByLogOccupancy.setOnAction((e)->root.setColorStyle(KStarTreeNode.ColorStyle.logOccupancy));
        MenuItem toggleCenter = new MenuItem("Show/Hide white center");
        toggleCenter.setOnAction(e->{
            root.toggleCenter();
        });
        help.getItems().add(helpDevShortCut3);
        file.getItems().add(loadTree);
        options.getItems().addAll(setvisibleLevels,toggleCenter, colorByEnergy, colorByOccupancy, colorByLogOccupancy);
        MenuBar menuBar = new MenuBar();
        menuBar.getMenus().addAll(file, options, help);
        Button button = new Button();
        button.setText("Click me!");
        return menuBar;
    }

    private void devShortCut3() {
        loadTreeFromFile(new File("Complex2XXMCATS.ordered.0.01.txt"));
    }

    private void devShortCut2() {
        loadTreeFromFile(new File("Protein2XXMContinuousBounds.txt"));
    }

    private void devShortCut() {
        loadTreeFromFile(new File("Complex2XXMContinuousBounds.txt"));
    }

    private void loadTreeFromFile(File selectedFile) {
        ringNode = new Pane();
        System.out.println("Parsing "+selectedFile);
        rootGroup = new Group();
        Group ringGroup = new Group();
        ringGroup.getChildren().add(dot);
        Group textGroup = new Group();
        Group g = rootGroup;
        Pane centerPane = new Pane();
        root = KStarTreeNode.parseTree(selectedFile, true);
        /*
        int level = 5;
        System.out.println("Enthalpy:"+root.computeEnthalpy(level));
        System.out.println("Entropy:"+root.computeEntropy(level));
        System.out.println("Num States at level "+level+":"+root.numStatesAtLevel(level));
        */
        rootGroup.getChildren().addAll(ringGroup, textGroup);
        root.setGroup(ringGroup);
        root.preprocess();
        root.render(g);
        root.setTextRoot(textGroup);
        root.autoExpand(0.001, 1);//,5);
        resize();
        //root.pieChart(1, 3,6);
        root.showRoot();
        centerPane.getChildren().addAll(g);
        triroot.setCenter(centerPane);
        centerPane.setOnScroll((event)-> {
            if (event.getDeltaY() == 0) {
                return;
            }

            double scaleFactor = (event.getDeltaY() > 0) ? SCALE_DELTA : 1 / SCALE_DELTA;
            double mouseX = event.getX();
            double mouseY = event.getY();
            Point2D mouseXY = new Point2D(mouseX, mouseY);
            Point2D mouseLocal = ringGroup.sceneToLocal(mouseXY);
            Point2D ringScene = ringGroup.localToScene(mouseXY);
            dot.setTranslateX(mouseLocal.getX());
            dot.setTranslateY(mouseLocal.getY());
            double ringWidth = ringGroup.getBoundsInLocal().getWidth();
            double ringHeight= ringGroup.getBoundsInLocal().getHeight();
            double ringCenterX = ringGroup.getTranslateX()-ringWidth/2;
            double ringCenterY = ringGroup.getTranslateY()-ringHeight/2;
            System.out.println("Mouse: ("+mouseX+","+mouseY+")");
            System.out.println("Mouse local: ("+mouseLocal.getX()+","+mouseLocal.getY()+")");

            System.out.println("Ring: ("+ringCenterX+","+ringCenterY+")");
            double distFromRingCenterX = ringCenterX-mouseLocal.getX();
            double distFromRingCenterY = ringCenterY-mouseLocal.getY();
            double deltaX = (1-scaleFactor)*distFromRingCenterX;
            double deltaY = (1-scaleFactor)*distFromRingCenterY;

            ringGroup.setScaleX(ringGroup.getScaleX() * scaleFactor);
            ringGroup.setScaleY(ringGroup.getScaleY() * scaleFactor);
            Point2D movedMouseScene = ringGroup.localToScene(mouseLocal);
            System.out.println("Moved mouse Scene: ("+movedMouseScene.getX()+","+movedMouseScene.getY()+")");
            ringGroup.setTranslateX(ringGroup.getTranslateX() + mouseX - movedMouseScene.getX());
            ringGroup.setTranslateY(ringGroup.getTranslateY() + mouseY - movedMouseScene.getY());
            //resize();
        });
        centerPane.setOnMousePressed((event)-> {
            mouseDownX = event.getX();
            mouseDownY = event.getY();
            Point2D mouseXY = new Point2D(mouseDownX, mouseDownY);
            Point2D mouseLocal = ringGroup.sceneToLocal(mouseXY);
            Point2D ringScene = ringGroup.localToScene(mouseXY);

            ringX = ringGroup.getTranslateX();
            ringY = ringGroup.getTranslateY();
        });
        centerPane.setOnMouseDragged((event)-> {
            double x = event.getX();
            double y = event.getY();
            ringGroup.setTranslateX(ringX+(x-mouseDownX));
            ringGroup.setTranslateY(ringY+(y-mouseDownY));

        });
        triroot.setTop(getMenuBar(primaryStage));
        //triroot.widthProperty().addListener(o-> resize());
        //triroot.heightProperty().addListener(o-> resize());
    }

    private void resize() {
        double width  = triroot.getWidth() - triroot.getInsets().getLeft() - triroot.getInsets().getRight();
        double height = triroot.getHeight() - triroot.getInsets().getTop() - triroot.getInsets().getBottom();

        triroot.setMaxSize(width, height);
        triroot.setPrefSize(width, height);
        rootGroup.setTranslateX(width/2);//-rootGroup.getBoundsInLocal().getWidth()/2);
        rootGroup.setTranslateY(height/2);//-rootGroup.getBoundsInLocal().getHeight()/2);

    }




    public static void main(String[] args) {
        launch(args);
    }
}
