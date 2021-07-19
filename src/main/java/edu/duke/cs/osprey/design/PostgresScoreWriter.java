package edu.duke.cs.osprey.design;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.SequenceAnalyzer;
import edu.duke.cs.osprey.tools.MathTools;
import org.sql2o.Sql2o;
import software.amazon.awssdk.core.sync.RequestBody;
import software.amazon.awssdk.regions.Region;
import software.amazon.awssdk.regions.RegionMetadata;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.S3ClientBuilder;
import software.amazon.awssdk.services.s3.model.PutObjectRequest;
import software.amazon.awssdk.services.s3.model.StorageClass;

import java.util.List;

public class PostgresScoreWriter implements KStarScoreWriter {

    private final String designName;
    private final List<String> inputs;
    private final Sql2o sql2o;
    private final int numConfsToSave;
    private final String s3BucketName;
    private final S3Client s3Client;
    int designId;

    public PostgresScoreWriter(PostgresConnectionInfo connectionInfo, S3Settings s3Settings,
                               String designName, List<String> inputs, int numConfsToSave) {
        this.designName = designName;
        this.inputs = inputs;
        this.numConfsToSave = numConfsToSave;
        this.s3Client = S3Client.builder().region(Region.of(s3Settings.region)).build();
        this.s3BucketName = s3Settings.bucketName;
        this.sql2o = new Sql2o(connectionInfo.connectionString, connectionInfo.username, connectionInfo.password);
    }

    @Override
    public void writeHeader() {
        var insertSql = "INSERT INTO designs (name, type) VALUES (:name, 'affinity')";

        try (var con = sql2o.open()) {
            designId = con.createQuery(insertSql, true)
                    .addParameter("name", designName)
                    .executeUpdate()
                    .getKey(Integer.class);
        }

        var inputSql = "INSERT INTO args (design, value) VALUES (:design, :value)";

        try (var con = sql2o.open()) {
            for (var input : inputs) {
                con.createQuery(inputSql, true)
                        .addParameter("design", designId)
                        .addParameter("value", input)
                        .executeUpdate();
            }
        }
    }

    @Override
    public void writeScore(ScoreInfo info) {
        var insertSql =
                "INSERT INTO affinities " +
                        "(design, is_wt, kstar_lower, kstar_upper, " +
                        "protein_lower, protein_upper, protein_confs_enumerated, protein_epsilon, " +
                        "ligand_lower, ligand_upper, ligand_confs_enumerated, ligand_epsilon, " +
                        "complex_lower, complex_upper, complex_confs_enumerated, complex_epsilon, " +
                        "variances) " +
                        "VALUES " +
                        "(:design, :is_wt, :kstar_lower, :kstar_upper, " +
                        ":protein_lower, :protein_upper, :protein_confs_enumerated, :protein_epsilon, " +
                        ":ligand_lower, :ligand_upper, :ligand_confs_enumerated, :ligand_epsilon, " +
                        ":complex_lower, :complex_upper, :complex_confs_enumerated, :complex_epsilon, " +
                        ":variances) ";

        var design = designId;
        var is_wt = info.sequenceNumber == 0;
        var kstar_lower = MathTools.log10(info.kstarScore.lowerBound);
        var kstar_upper = MathTools.log10(info.kstarScore.upperBound);

        var protein_lower = MathTools.log10(info.kstarScore.protein.values.calcLowerBound());
        var protein_upper = MathTools.log10(info.kstarScore.protein.values.calcUpperBound());
        var protein_confs_enumerated = info.kstarScore.protein.numConfs;
        var protein_epsilon = info.kstarScore.protein.values.getEffectiveEpsilon();

        var ligand_lower = MathTools.log10(info.kstarScore.ligand.values.calcLowerBound());
        var ligand_upper = MathTools.log10(info.kstarScore.ligand.values.calcUpperBound());
        var ligand_confs_enumerated = info.kstarScore.ligand.numConfs;
        var ligand_epsilon = info.kstarScore.ligand.values.getEffectiveEpsilon();

        var complex_lower = MathTools.log10(info.kstarScore.complex.values.calcLowerBound());
        var complex_upper = MathTools.log10(info.kstarScore.complex.values.calcUpperBound());
        var complex_confs_enumerated = info.kstarScore.complex.numConfs;
        var complex_epsilon = info.kstarScore.complex.values.getEffectiveEpsilon();

        var variances = info.sequence.toString(Sequence.Renderer.AssignmentMutations, info.sequence.calcCellSize() + 1).trim();

        int affinityId;
        try (var con = sql2o.open()) {
            affinityId =
                con.createQuery(insertSql, true)
                    .addParameter("design", design)
                    .addParameter("is_wt", is_wt)
                    .addParameter("kstar_lower", kstar_lower)
                    .addParameter("kstar_upper", kstar_upper)
                    .addParameter("protein_lower", protein_lower)
                    .addParameter("protein_upper", protein_upper)
                    .addParameter("protein_confs_enumerated", protein_confs_enumerated)
                    .addParameter("protein_epsilon", protein_epsilon)
                    .addParameter("ligand_lower", ligand_lower)
                    .addParameter("ligand_upper", ligand_upper)
                    .addParameter("ligand_confs_enumerated", ligand_confs_enumerated)
                    .addParameter("ligand_epsilon", ligand_epsilon)
                    .addParameter("complex_lower", complex_lower)
                    .addParameter("complex_upper", complex_upper)
                    .addParameter("complex_confs_enumerated", complex_confs_enumerated)
                    .addParameter("complex_epsilon", complex_epsilon)
                    .addParameter("variances", variances)
                    .executeUpdate()
                    .getKey(Integer.class);
        }

        if (numConfsToSave > 0) {
            var analysis = new SequenceAnalyzer(info.kstar).analyze(info.sequence, numConfsToSave);
            var pdb = analysis.ensemble.writePdbString(String.format("Top %d confs for sequence", numConfsToSave));

            var pdbSql = "INSERT INTO affinity_structures (affinity, structure) VALUES (:affinity, :structure)";

            int structId;
            try (var con = sql2o.open()) {
                structId = con.createQuery(pdbSql, true)
                        .addParameter("affinity", affinityId)
                        .addParameter("structure", "")
                        .executeUpdate()
                        .getKey(Integer.class);
            }

            var s3Key = String.format("%s-%s.pdb", designName, structId);
            var request = PutObjectRequest.builder()
                    .bucket(s3BucketName)
                    .storageClass(StorageClass.INTELLIGENT_TIERING)
                    .key(s3Key)
                    .build();

            var response = s3Client.putObject(request, RequestBody.fromString(pdb));

            var structureUpdateSql = "UPDATE affinity_structures set structure = :s3Key where id = :structId";
            try (var con = sql2o.open()) {
                con.createQuery(structureUpdateSql, false)
                        .addParameter("s3Key", /*s3://{bucket}/{key}*/String.format("s3://%s/%s", s3BucketName, s3Key))
                        .addParameter("structId", structId)
                        .executeUpdate();
            }
        }
    }
}
