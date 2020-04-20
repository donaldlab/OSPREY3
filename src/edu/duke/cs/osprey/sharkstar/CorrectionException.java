package edu.duke.cs.osprey.sharkstar;

public class CorrectionException extends RuntimeException{
    public MultiSequenceSHARKStarNode offendingNode;
    public double HOTCorrection;

    public CorrectionException(){
            super();
        }

    public CorrectionException(String message) {
        super(message);
    }

    public CorrectionException(String message, Throwable cause) {
        super(message, cause);
    }
    public CorrectionException(Throwable cause) {
        super(cause);
    }

    public CorrectionException(String message, MultiSequenceSHARKStarNode offendingNode, double HOTCorrection){
        super(message.concat(
                String.format("\nNode: %s, Correction: %.3f",
                        offendingNode.getConfSearchNode().confToString(), HOTCorrection)
        ));
        this.offendingNode = offendingNode;
        this.HOTCorrection = HOTCorrection;

    }

    protected CorrectionException(String message, Throwable cause,
                                  boolean enableSuppression,
                                  boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }

    public void setOffendingNode(MultiSequenceSHARKStarNode node){
        this.offendingNode = node;
    }
    public void setHOTCorrection(double corr){
        this.HOTCorrection = corr;
    }

}
