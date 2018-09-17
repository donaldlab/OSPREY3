package edu.duke.cs.osprey.markstar.visualizer;

import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.stage.FileChooser;
import javafx.stage.Stage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

public class Visualizer extends Application {

    private static final double SCALE_DELTA = 1.1;
    Stage primaryStage;
    BorderPane triroot;
    Pane ringNode;
    KStarTreeNode root;
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
        help.getItems().add(helpDevShortCut3);
        file.getItems().add(loadTree);
        MenuBar menuBar = new MenuBar();
        menuBar.getMenus().addAll(file, options, help);
        Button button = new Button();
        button.setText("Click me!");
        triroot.setTop(menuBar);
        primaryStage.setScene(test);
        primaryStage.show();

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
        Group textGroup = new Group();
        Group g = rootGroup;
        Pane centerPane = new Pane();
        root = KStarTreeNode.parseTree(selectedFile, true);
        int level = 5;
        System.out.println("Enthalpy:"+root.computeEnthalpy(level));
        System.out.println("Entropy:"+root.computeEntropy(level));
        System.out.println("Num States at level "+level+":"+root.numStatesAtLevel(level));
        rootGroup.getChildren().addAll(ringGroup, textGroup);
        root.setGroup(ringGroup);
        root.preprocess();
        root.render(g);
        root.setTextRoot(textGroup);
        root.autoExpand(0.001);//,5);
        resize();
        //root.pieChart(1);
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
            double ringWidth = ringGroup.getBoundsInLocal().getWidth();
            double ringHeight= ringGroup.getBoundsInLocal().getHeight();
            double ringCenterX = ringGroup.getTranslateX()+ringWidth/2;
            double ringCenterY = ringGroup.getTranslateX()+ringHeight/2;
            double distFromRingCenterX = ringCenterX-mouseX;
            double distFromRingCenterY = ringCenterY-mouseY;
            double deltaX = (scaleFactor-1)*distFromRingCenterX;
            double deltaY = (scaleFactor-1)*distFromRingCenterY;

            ringGroup.setScaleX(ringGroup.getScaleX() * scaleFactor);
            ringGroup.setScaleY(ringGroup.getScaleY() * scaleFactor);
            ringGroup.setTranslateX(-deltaX);
            ringGroup.setTranslateX(-deltaY);
            resize();
        });
        centerPane.setOnMousePressed((event)-> {
            mouseDownX = event.getX();
            mouseDownY = event.getY();
            ringX = ringGroup.getTranslateX();
            ringY = ringGroup.getTranslateY();
        });
        centerPane.setOnMouseDragged((event)-> {
            double x = event.getX();
            double y = event.getY();
            ringGroup.setTranslateX(ringX+(x-mouseDownX));
            ringGroup.setTranslateY(ringY+(y-mouseDownY));

        });
        final Menu file = new Menu("File");
        final Menu options = new Menu("Options");
        final Menu help = new Menu("Help");
        MenuItem loadTree = new MenuItem("Load tree!");
        loadTree.setOnAction(e -> {
            FileChooser fc = new FileChooser();
            File selectedFile2 = fc.showOpenDialog(primaryStage);
            loadTreeFromFile(selectedFile2);
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
        help.getItems().add(helpDevShortCut3);
        file.getItems().add(loadTree);
        MenuBar menuBar = new MenuBar();
        menuBar.getMenus().addAll(file, options, help);
        triroot.setTop(menuBar);
        triroot.widthProperty().addListener(o-> resize());
        triroot.heightProperty().addListener(o-> resize());
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
