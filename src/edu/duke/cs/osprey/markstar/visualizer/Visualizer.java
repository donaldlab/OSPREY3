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

    Stage primaryStage;
    BorderPane triroot;
    Pane ringNode;
    KStarTreeNode root;
    Group rootGroup;

    @Override
    public void start(Stage primaryStage) throws Exception{
        this.primaryStage = primaryStage;
        primaryStage.setTitle("MARK* Tree Analyzer (v0.2)");
        BorderPane root = new BorderPane();
        Scene test = new Scene(root, 300, 275);
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
        root.setTop(menuBar);
        primaryStage.setScene(test);
        primaryStage.show();

    }

    private void devShortCut3() {
        loadTreeFromFile(new File("Complex2XXMCATS.ordered.0.01.txt"));
    }

    private void devShortCut2() {
        loadTreeFromFile(new File("ProteinConfTreeBounds.txt"));
    }

    private void devShortCut() {
        loadTreeFromFile(new File("ComplexConfTreeBounds.txt"));
    }

    private void loadTreeFromFile(File selectedFile) {
        triroot = new BorderPane();
        triroot.setBackground(new Background(new BackgroundFill(Color.color(1,1,1), CornerRadii.EMPTY, Insets.EMPTY)));
        ringNode = new Pane();
        Scene tri = new Scene(triroot, primaryStage.getWidth(), primaryStage.getHeight());
        System.out.println("Parsing "+selectedFile);
        rootGroup = new Group();
        Group g = rootGroup;
        ScrollPane centerPane = new ScrollPane();
        root = KStarTreeNode.parseTree(selectedFile, true);
        root.setGroup(g);
        root.preprocess();
        root.render(g);
        root.autoExpand(0.01);
        resize();
        root.showRoot();
        centerPane.setContent(g);
        triroot.setCenter(centerPane);
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
        primaryStage.setScene(tri);
        triroot.widthProperty().addListener(o-> resize());
        triroot.heightProperty().addListener(o-> resize());
    }

    private void resize() {
        double width  = triroot.getWidth() - triroot.getInsets().getLeft() - triroot.getInsets().getRight();
        double height = triroot.getHeight() - triroot.getInsets().getTop() - triroot.getInsets().getBottom();

        triroot.setMaxSize(width, height);
        triroot.setPrefSize(width, height);
        rootGroup.setTranslateX(width/2-rootGroup.getBoundsInLocal().getWidth()/2);
        rootGroup.setTranslateY(height/2-rootGroup.getBoundsInLocal().getHeight()/2);

    }




    public static void main(String[] args) {
        launch(args);
    }
}
