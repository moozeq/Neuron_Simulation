#pragma once
#include "Utilities.h"

const float YAW = -90.0f;
const float PITCH = 0.0f;
const float SPEED = 0.5f;
const float SENSITIVTY = 0.07f;
const float ZOOM = 45.0f;

namespace cam {
	enum CameraMovement {
		FORWARD,
		BACKWARD,
		LEFT,
		RIGHT
	};
}

class Camera
{
	friend class Simulation;
	glm::vec3 position;
	glm::vec3 front;
	glm::vec3 up;
	glm::vec3 right;
	glm::vec3 worldUp;

	float yaw;
	float pitch;

	float movementSpeed;
	float mouseSensitivity;
	float zoom;

	float lastX;
	float lastY;

	void updateCameraVectors() {
		glm::vec3 front;
		front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
		front.y = sin(glm::radians(pitch));
		front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
		this->front = glm::normalize(front);
		right = glm::normalize(glm::cross(this->front, worldUp));
		up = glm::normalize(glm::cross(right, this->front));
	}

	Camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f), float yaw = YAW, float pitch = PITCH) :
		front(glm::vec3(0.0f, 0.0f, -1.0f)), movementSpeed(SPEED), mouseSensitivity(SENSITIVTY), zoom(ZOOM)
	{
		this->position = position;
		this->worldUp = up;
		this->yaw = yaw;
		this->pitch = pitch;
		updateCameraVectors();
	}

	Camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) :
		front(glm::vec3(0.0f, 0.0f, -1.0f)), movementSpeed(SPEED), mouseSensitivity(SENSITIVTY), zoom(ZOOM)
	{
		this->position = glm::vec3(posX, posY, posZ);
		this->worldUp = glm::vec3(upX, upY, upZ);
		this->yaw = yaw;
		this->pitch = pitch;
		updateCameraVectors();
	}

	void processKeyboard(cam::CameraMovement direction, GLfloat deltaTime) {
		float velocity = movementSpeed * deltaTime;
		if (direction == cam::FORWARD)
			position += front * velocity;
		if (direction == cam::BACKWARD)
			position -= front * velocity;
		if (direction == cam::LEFT)
			position -= right * velocity;
		if (direction == cam::RIGHT)
			position += right * velocity;
	}

	void processMouseMovement(GLfloat xoffset, GLfloat yoffset, GLboolean constrainPitch = GL_TRUE) {
		xoffset *= mouseSensitivity;
		yoffset *= mouseSensitivity;

		yaw += xoffset;
		pitch += yoffset;

		if (constrainPitch)
		{
			if (pitch > 89.0f)
				pitch = 89.0f;
			if (pitch < -89.0f)
				pitch = -89.0f;
		}
		updateCameraVectors();
	}

	void processMouseScroll(GLfloat yoffset) {
		if (zoom >= 1.0f && zoom <= 45.0f)
			zoom -= yoffset;
		if (zoom <= 1.0f)
			zoom = 1.0f;
		if (zoom >= 45.0f)
			zoom = 45.0f;
	}

	glm::mat4 GetViewMatrix()
	{
		return glm::lookAt(position, position + front, up);
	}
};