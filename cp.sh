#!/bin/bash
mount /dev/sdf1 ../boot
cp arch/arm64/boot/Image ../boot
umount /dev/sdf1
