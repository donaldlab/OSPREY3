package edu.duke.cs.osprey.gui.view

import cuchaz.kludge.vulkan.ColorRGBA
import cuchaz.kludge.vulkan.putColor4Bytes
import edu.duke.cs.osprey.molscope.render.CylinderRenderable
import edu.duke.cs.osprey.molscope.render.RenderEffect
import edu.duke.cs.osprey.molscope.render.SphereRenderable
import edu.duke.cs.osprey.molscope.render.put
import edu.duke.cs.osprey.molscope.view.Color
import edu.duke.cs.osprey.molscope.view.ColorPalette
import edu.duke.cs.osprey.molscope.view.ColorsMode
import edu.duke.cs.osprey.molscope.view.RenderView
import org.joml.Vector3d
import java.nio.ByteBuffer


class ProbeView : RenderView {

	data class Group(
		val id: String,
		val dots: Map<String,List<Vector3d>>,
		val vectors: Map<String,List<Pair<Vector3d,Vector3d>>>
	)

	override var isVisible = true

	private var sequence = 0

	inner class Visibility {

		private val flags = HashMap<String,Boolean>()

		operator fun get(id: String) =
			flags.getOrPut(id) { true }

		operator fun set(id: String, value: Boolean) {
			flags[id] = value
			sequence += 1
		}
	}
	val visibility = Visibility()

	var groups: Map<String,Group>? = null
		set(value) {
			field = value
			sequence += 1
		}

	var dotRadius: Double = 0.04
		set(value) {
			field = value
			sequence += 1
		}

	var vectorRadius: Double = 0.02
		set(value) {
			field = value
			sequence += 1
		}

	private val visibleGroups get() =
		groups?.values
			?.filter { visibility[it.id] }
			?: emptyList()

	// render the "dots" as spheres
	override val spheres = object : SphereRenderable {

        override val numVertices
            get() = visibleGroups
                .sumOf { group -> group.dots.values.sumOf { it.size * 4 } }

        override val verticesSequence get() = sequence

        override fun fillVertexBuffer(buf: ByteBuffer, colorsMode: ColorsMode) {

            val radius = dotRadius.toFloat()

            for (group in visibleGroups) {
                for ((colorName, dots) in group.dots) {

                    val color = getColor(colorName)[colorsMode]

                    for (dot in dots) {

                        // copy the vertices 4 times
                        for (i in 0 until 4) {

                            // downgrade pos to floats for rendering
                            buf.putFloat(dot.x.toFloat())
                            buf.putFloat(dot.y.toFloat())
                            buf.putFloat(dot.z.toFloat())

                            buf.putFloat(radius)
                            buf.putColor4Bytes(color)

                            // no effects needed (yet?)
                            buf.putInt(0)
                            buf.put(null as RenderEffect?)
                        }
                    }
                }
            }
        }
    }

	// render the "vectors" as cylinders
	override val cylinders = object : CylinderRenderable {

        override val numVertices
            get() = visibleGroups
                .sumOf { group -> group.vectors.values.sumOf { it.size * 4 } }

        override val verticesSequence get() = sequence

        override fun fillVertexBuffer(buf: ByteBuffer, colorsMode: ColorsMode) {

            val radius = vectorRadius.toFloat()

            for (group in visibleGroups) {
                for ((colorName, vectors) in group.vectors) {
                    val color = getColor(colorName)[colorsMode]

                    for ((a, b) in vectors) {
                        val points = listOf(a, b)

                        // write each cylinder 4 times
                        for (i in 0 until 4) {

                            for (p in points) {
                                // downgrade pos to floats for rendering
                                buf.putFloat(p.x.toFloat())
                                buf.putFloat(p.y.toFloat())
                                buf.putFloat(p.z.toFloat())
                            }

                            for (p in points) {
                                buf.putFloat(radius)
                            }

                            for (p in points) {
                                buf.putColor4Bytes(color)
                            }

                            for (p in points) {
                                buf.putInt(0)
                            }

                            // no effects needed (yet?)
                            for (p in points) {
                                buf.put(null as RenderEffect?)
                            }
                        }
                    }
                }
            }
        }
    }

	companion object {

		fun getColor(name: String): Color =
			colors[name] ?: ColorPalette.lightGrey

		// kinemage colors
		val colors: Map<String,Color> = HashMap<String,Color>().apply {

			fun hsv(name: String, hue: Int, darkSat: Int, lightSat: Int, darkVal: Int, lightVal: Int) {
				this[name] = Color(
					name,
					dark = java.awt.Color.getHSBColor(
						hue.toFloat()/360f,
						darkSat.toFloat()/100f,
						darkVal.toFloat()/100f
					).run {
						ColorRGBA.Int(red, green, blue)
					},
					light = java.awt.Color.getHSBColor(
						hue.toFloat()/360f,
						lightSat.toFloat()/100f,
						lightVal.toFloat()/100f
					).run {
						ColorRGBA.Int(red, green, blue)
					}
				)
			}

			hsv("default",   0,      0,      0,      100,    0)
			hsv("red",       0,      100,    100,    100,    80)
			hsv("orange",    20,     100,    100,    100,    90)
			hsv("rust",      20,     100,    100,    100,    90)
			hsv("gold",      40,     100,    100,    100,    90)
			hsv("yellow",    60,     100,    100,    100,    90)
			hsv("lime",      80,     100,    100,    100,    85)
			hsv("green",     120,    80,     90,     100,    75)
			hsv("sea",       150,    100,    100,    100,    85)
			hsv("seagreen",  150,    100,    100,    100,    85)
			hsv("cyan",      180,    100,    85,     85,     80)
			hsv("sky",       210,    75,     80,     95,     90)
			hsv("skyblue",   210,    75,     80,     95,     90)
			hsv("blue",      240,    70,     80,     100,    100)
			hsv("purple",    275,    75,     100,    100,    85)
			hsv("magenta",   300,    95,     100,    100,    90)
			hsv("hotpink",   335,    100,    100,    100,    90)
			hsv("pink",      350,    55,     75,     100,    90)
			hsv("peach",     25,     75,     75,     100,    90)
			hsv("lilac",     275,    55,     75,     100,    80)
			hsv("pinktint",  340,    30,     100,    100,    55)
			hsv("peachtint", 25,     50,     100,    100,    60)
			hsv("yellowtint",60,     50,     100,    100,    75)
			hsv("paleyellow",60,     50,     100,    100,    75)
			hsv("greentint", 135,    40,     100,    100,    35)
			hsv("bluetint",  220,    40,     100,    100,    50)
			hsv("lilactint", 275,    35,     100,    100,    45)
			hsv("white",     0,      0,      0,      100,    0)
			hsv("gray",      0,      0,      0,      50,     40)
			hsv("grey",      0,      0,      0,      50,     40)
			hsv("brown",     20,     45,     45,     75,     55)
		}
	}
}