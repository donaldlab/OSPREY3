package edu.duke.cs.osprey.gui.tools


/**
 * An array list that implements a map keyed by integers.
 * Very efficient map implementation for get/put when using consecutive keys in [0,N).
 */
class ArrayMap<T:Any> : MutableMap<Int,T> {

	private val list = ArrayList<T?>()

	override var size: Int = 0
		private set

	override fun containsKey(key: Int): Boolean =
		this[key] != null

	override fun containsValue(value: T): Boolean =
		list.contains(value)

	override fun get(key: Int): T? =
		if (key >= 0 && key < list.size) {
			list[key]
		} else {
			null
		}

	override fun isEmpty(): Boolean =
		size == 0

	private inner class KeyIterator : Iterator<Int> {

		private var i: Int? = findNext(-1)

		fun findNext(i: Int): Int? {
			for (j in i+1 until list.size) {
				if (list[j] != null) {
					return j
				}
			}
			return null
		}

		override fun hasNext(): Boolean =
			i != null

		override fun next(): Int {
			val i = i ?: throw NoSuchElementException()
			this.i = findNext(i)
			return i
		}
	}

	private inner class Entry(
		override val key: Int
	) : MutableMap.MutableEntry<Int, T> {

		private val map = this@ArrayMap

		override val value: T get() = map[key]!!

		override fun setValue(newValue: T): T =
			map.put(key, newValue)!!
	}

	override val entries: MutableSet<MutableMap.MutableEntry<Int, T>> by lazy {
		object : MutableSet<MutableMap.MutableEntry<Int, T>> {

			val map get() = this@ArrayMap

			override fun add(element: MutableMap.MutableEntry<Int, T>): Boolean {
				map[element.key] = element.value
				return true
			}

			override fun addAll(elements: Collection<MutableMap.MutableEntry<Int, T>>): Boolean {
				for (element in elements) {
					add(element)
				}
				return true
			}

			override fun clear() {
				map.clear()
			}

			override fun iterator(): MutableIterator<MutableMap.MutableEntry<Int, T>> =
				object : MutableIterator<ArrayMap<T>.Entry> {

					val keys = KeyIterator()
					var key: Int? = null

					override fun hasNext() =
						keys.hasNext()

					override fun next() =
						Entry(keys.next())
							.also { key = it.key }

					override fun remove() {
						key?.let { map.remove(it) }
							?: throw NoSuchElementException()
					}
				}

			override fun remove(element: MutableMap.MutableEntry<Int, T>): Boolean =
				map.remove(element.key) != null

			override fun removeAll(elements: Collection<MutableMap.MutableEntry<Int, T>>) =
				elements
					.map { remove(it) }
					.any { it }

			override fun retainAll(elements: Collection<MutableMap.MutableEntry<Int, T>>) =
				throw UnsupportedOperationException()

			override val size
				get() = map.size

			override fun contains(element: MutableMap.MutableEntry<Int, T>) =
				map.containsKey(element.key)

			override fun containsAll(elements: Collection<MutableMap.MutableEntry<Int, T>>) =
				elements.all { contains(it) }

			override fun isEmpty() =
				map.isEmpty()
		}
	}

	override val keys: MutableSet<Int> by lazy {
		object : MutableSet<Int> {

			private val map get() = this@ArrayMap

			override fun add(element: Int) =
				throw UnsupportedOperationException()

			override fun addAll(elements: Collection<Int>) =
				throw UnsupportedOperationException()

			override fun clear() =
				map.clear()

			override fun iterator(): MutableIterator<Int> =
				object : MutableIterator<Int> {

					val keys = KeyIterator()
					var key: Int? = null

					override fun hasNext() =
						keys.hasNext()

					override fun next() =
						keys.next()
							.also { key = it }

					override fun remove() {
						key?.let { map.remove(it) }
							?: throw NoSuchElementException()
					}
				}

			override fun remove(element: Int) =
				map.remove(element) != null

			override fun removeAll(elements: Collection<Int>) =
				elements
					.map { remove(it) }
					.any { it }

			override fun retainAll(elements: Collection<Int>) =
				throw UnsupportedOperationException()

			override val size
				get() = map.size

			override fun contains(element: Int) =
				map.containsKey(element)

			override fun containsAll(elements: Collection<Int>) =
				elements.all { contains(it) }

			override fun isEmpty() =
				map.isEmpty()
		}
	}

	override val values: MutableCollection<T> by lazy {
		object : MutableCollection<T> {

			val map get() = this@ArrayMap

			override val size =
				map.size

			override fun contains(element: T) =
				map.containsValue(element)

			override fun containsAll(elements: Collection<T>) =
				elements.all { contains(it) }

			override fun isEmpty() =
				map.isEmpty()

			override fun add(element: T) =
				throw UnsupportedOperationException()

			override fun addAll(elements: Collection<T>) =
				throw UnsupportedOperationException()

			override fun clear() =
				map.clear()

			override fun iterator(): MutableIterator<T> =
				object : MutableIterator<T> {

					val keys = KeyIterator()
					var key: Int? = null

					override fun hasNext() =
						keys.hasNext()

					override fun next() =
						keys.next().let { i ->
							key = i
							list[i]!!
						}

					override fun remove() {
						key?.let { map.remove(it) }
							?: throw NoSuchElementException()
					}
				}

			override fun remove(element: T) =
				throw UnsupportedOperationException()

			override fun removeAll(elements: Collection<T>) =
				throw UnsupportedOperationException()

			override fun retainAll(elements: Collection<T>) =
				throw UnsupportedOperationException()
		}
	}

	override fun clear() =
		list.clear()

	override fun put(key: Int, value: T): T? {
		if (key < 0) {
			throw IllegalArgumentException("keys cann't be less than zero")
		}
		while (key >= list.size) {
			list.add(null)
		}
		val old = list[key]
		list[key] = value
		if (old == null) {
			size += 1
		}
		return old
	}

	override fun putAll(from: Map<out Int, T>) {
		for ((key, value) in from) {
			put(key, value)
		}
	}

	override fun remove(key: Int): T? {
		val old = this[key]
		if (old != null) {
			list[key] = null
			size -= 1
		}
		return old
	}
}
