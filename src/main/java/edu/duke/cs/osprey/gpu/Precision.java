package edu.duke.cs.osprey.gpu;

/**
 * Determines the level of floating-point precision for calculations.
 */
public enum Precision {

    Float32(4) {

        @Override
        public Object fromDouble(double val) {
            return (float)val;
        }

        @Override
        public double toDouble(Object val) {
            return (double)(Float)val;
        }
    },

    Float64(8) {

        @Override
        public Object fromDouble(double val) {
            return val;
        }

        @Override
        public double toDouble(Object val) {
            return (Double)val;
        }
    };

    public final int bytes;

    Precision(int bytes) {
        this.bytes = bytes;
    }

    public abstract Object fromDouble(double val);
    public abstract double toDouble(Object val);

    public double cast(double val) {
        return toDouble(fromDouble(val));
    }

    public <T> T map(T f32, T f64) {
        return switch (this) {
            case Float32 -> f32;
            case Float64 -> f64;
        };
    }
}
