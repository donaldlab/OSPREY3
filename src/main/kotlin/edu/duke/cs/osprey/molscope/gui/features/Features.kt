package edu.duke.cs.osprey.molscope.gui.features


class Features<T:HasFeatureId> {

	inner class Menu(val name: String) {

		val id = FeatureId(name)

		val features: List<T> get() = _features
		private val _features = ArrayList<T>()

		fun find(id: FeatureId) = _features.find { it.id == id }
		fun contains(id: FeatureId) = find(id) != null

		fun add(feature: T) {

			// check for duplicates
			if (contains(feature.id)) {
				throw IllegalArgumentException("feature already exists in this menu")
			}

			_features.add(feature)
		}
	}

	val menus: List<Menu> get() = _menus
	private val _menus = ArrayList<Menu>()

	fun menu(name: String): Menu {
		val id = FeatureId(name)
		return _menus.find { it.id == id }
			?: Menu(name).apply {
				_menus.add(this)
			}
	}

	fun <R> menu(name: String, block: Menu.() -> R): R = menu(name).block()
}
