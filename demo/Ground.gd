extends Node

func _on_GH_value_changed(value):
	self.translation = Vector3(0, value, 0);
