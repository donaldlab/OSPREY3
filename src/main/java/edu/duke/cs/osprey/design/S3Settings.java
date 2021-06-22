package edu.duke.cs.osprey.design;

public class S3Settings {
    public String region;
    public String bucketName;

    public S3Settings(String region, String bucketName) {
        this.region = region;
        this.bucketName = bucketName;
    }
}
