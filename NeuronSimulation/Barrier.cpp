#include "Barrier.h"

Barrier::Barrier() :
	x0(0.0f), y0(0.0f), z0(0.0f)
{
	bilayerTexture = loadMipmapTexture(GL_TEXTURE0, "bilayer.png");
}

void Barrier::render() const
{
	glBindTexture(GL_TEXTURE_2D, bilayerTexture);
	glBindVertexArray(innerLayerVAO);
	glDrawElements(GL_TRIANGLES, innerLayerIndices.size() * 3, GL_UNSIGNED_INT, 0);
	glBindVertexArray(outerLayerVAO);
	glDrawElements(GL_TRIANGLES, outerLayerIndices.size() * 3, GL_UNSIGNED_INT, 0);
}
