package edu.duke.cs.osprey.sharkstar;

public class UpperBoundException extends RuntimeException{
    public MultiSequenceSHARKStarNode offendingNode;

    public UpperBoundException(){
        super();
    }
    public UpperBoundException(String message) {
        super(message);
    }

    public UpperBoundException(String message, Throwable cause) {
            super(message, cause);
        }
    public UpperBoundException(Throwable cause) {
        super(cause);
    }

    protected UpperBoundException(String message, Throwable cause,
                               boolean enableSuppression,
                               boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }

    public void setOffendingNode(MultiSequenceSHARKStarNode node){
        this.offendingNode = node;
    }

}
