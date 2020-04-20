package edu.duke.cs.osprey.sharkstar;

public class LowerBoundException extends RuntimeException{
    public MultiSequenceSHARKStarNode offendingNode;

    public LowerBoundException(){
        super();
    }
    public LowerBoundException(String message) {
        super(message);
    }

    public LowerBoundException(String message, Throwable cause) {
        super(message, cause);
    }
    public LowerBoundException(Throwable cause) {
        super(cause);
    }

    protected LowerBoundException(String message, Throwable cause,
                                  boolean enableSuppression,
                                  boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }

    public void setOffendingNode(MultiSequenceSHARKStarNode node){
        this.offendingNode = node;
    }

}
